/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ecfeci_sys.h>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <cal_leap_seconds.h>
#include <astro_eop_sys.h>

#include <sofa.h>
#include <sofam.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <memory>
#include <stdexcept>

/*
 * Local utility for converting a C double[3][3] to an Eigen Matrix3d.
 * Should hopefully be the only need for this utility in the codebase
 * unless other additional libraries, similar to SOFA, rely on C array
 * matrix representations.
 */
static Eigen::Matrix3d from3x3(double mtx[3][3]);

/*
 * Nomenclature (IAU typically uses "System" vs. "Frame", e.g., GCRS
 * vs.. GCRF.
 *   ITRF  International terrestrial reference frame (ECEF)
 *   TIRF  Terrestrial intermediate reference frame (PEF)
 *   CIRF  Celestial intermediate reference frame (non-equinox based TETE)
 *   GCRF  Geocentric  celestial reference frame (ECI, within ~20 mas of J2000)
 */

namespace {
    // Default buffer from jdStart and jdStop with a zero rate_days
  constexpr double offset_0 {0.25};
    // Percentage of rate_days to use as buffer for jdStart and jdStop
  constexpr double p_offset {0.25};
}

namespace eom {

EcfEciSys::EcfEciSys(const JulianDate& startTime,
                     const JulianDate& stopTime,
                     const Duration& dt,
                     const std::shared_ptr<eom::EopSys>& eopSys,
                     bool interpolate) :
                     jdStart {startTime}, jdStop {stopTime},
                     rate_days {dt.getDays()}, interpolate_bpnpm {interpolate}
{
    // Put a small buffer around start and stop to minimize logic
    // locating ECFECI data
  if (rate_days == 0.0) {
    double offset {offset_0};
    jdStart += -offset;
    jdStop  +=  offset;
  } else {
    double offset {p_offset*rate_days};
    jdStart += -offset;
      // Since backing off start time, need to push end time a
      // full message forward to not end up short
    jdStop  +=  rate_days;
  }

    // Variables used for intermediate calculations within loop.
  double x;                       // Celestial Intermediate Pole x coordinate
  double y;                       // Celestial Intermediate Pole y coordinate
  double s;                       // CIO locator, radians
  double cirf2gcrf[3][3];         // CIRF to GCRF
  double itrf2tirf[3][3];         // ITRF to TIRF
  double mod2j2000[3][3];         // MOD to J2000
    //
  LeapSeconds& ls = LeapSeconds::getInstance();
  auto jd = jdStart;
    // If only one message, compute for middle time and set rate
    // to be the full duration of the time period
  if (rate_days == 0.0) {
    rate_days = jdStop - jdStart;
    jd += rate_days/2.0;
  }
  while (jd <= jdStop) {
    ecf_eci f2i;
    meme_eci i2i;
    f2i.mjd2000 = jd.getMjd2000();
    i2i.mjd2000 = f2i.mjd2000;
    eop_record eop;
    auto jdTT = ls.utc2tt(jd);
      // Pole locations for BPN
    iauXys06a(jdTT.getJdHigh(), jdTT.getJdLow(), &x, &y, &s);
      // Get EOP data is available
    if (eopSys != nullptr) {
      eop = eopSys->getEop(jd);
      f2i.ut1mutc = phy_const::tu_per_sec*eop.ut1mutc;
      f2i.lod = phy_const::tu_per_sec*0.001*eop.lod;       // TU per msec
      x += utl_const::rad_per_mas*eop.dx;
      y += utl_const::rad_per_mas*eop.dy;
    }
      // Transpose of BPN, to BPN
    iauC2ixys(x, y, s, cirf2gcrf); 
    iauTr(cirf2gcrf, cirf2gcrf);
      // CIO locator used for polar motion transformation
    double sp {iauSp00(jdTT.getJdHigh(), jdTT.getJdLow())};
      // Polar motion, ECF to ECI after transpose
    iauPom00(utl_const::rad_per_arcsec*eop.xp,
             utl_const::rad_per_arcsec*eop.yp, sp, itrf2tirf);
    iauTr(itrf2tirf, itrf2tirf);
      // Convert to final form and insert
    Eigen::Matrix3d mbpn = from3x3(cirf2gcrf);
    Eigen::Matrix3d mpm = from3x3(itrf2tirf);
    Eigen::Quaterniond qbpn(mbpn);
    Eigen::Quaterniond qpm(mpm);
    f2i.pm = qpm;
    f2i.bpn = qbpn;
      // IAU 76 Precession to support legacy celestial algorithms
    iauPmat76(jdTT.getJdHigh(), jdTT.getJdLow(), mod2j2000);
    iauTr(mod2j2000, mod2j2000);
    Eigen::Matrix3d mp76 = from3x3(mod2j2000);
    Eigen::Quaterniond qp76(mp76);
    i2i.p76 = qp76;
      // Store records
    f2iData.push_back(f2i);
    memeData.push_back(i2i);

    jd += rate_days;
  }
    // Frame Bias - time independent
  {
    double gamb, phib, psib, eps;
      // J2000 to GCRF
    double rbt[3][3];
      // Bias-precession - first compute GCRF to J2000
    iauPfw06(DJM0, DJM00, &gamb, &phib, &psib, &eps);
    iauFw2m(gamb, phib, psib, eps, rbt);
      // Then the transpose for J2000 to GCRF
    iauTr(rbt, rbt);
    Eigen::Matrix3d mrbt = from3x3(rbt);
    Eigen::Quaterniond qrbt(mrbt);
    bt = qrbt;
  }
    // Size fixed after initialization
  nfi = static_cast<unsigned long>(f2iData.size());
}


ecf_eci EcfEciSys::getEcfEciData(const JulianDate& utc) const
{
    // Check for valid date
  double days {utc - jdStart};
  if (days < 0.0  ||  jdStop - utc < 0.0) {
    throw std::out_of_range("EcfEciSys::getEcfEciData() Time out of range");
  }

    // Always need first index
  unsigned long int ndx1 {static_cast<unsigned long int>(days/rate_days)};
  const ecf_eci& f2i1 = f2iData[ndx1];

    // Get second data set if interpolating - otherwise,
    // return data less than or equal to requested time
  if (nfi == 1UL) {
    return f2iData[0UL];
  } else if (interpolate_bpnpm) {
    double mjd2000 {utc.getMjd2000()};
    double dt_days {mjd2000 - f2i1.mjd2000};
    double dt {dt_days/rate_days};
    unsigned long int ndx2 {ndx1 + 1UL};
    const ecf_eci& f2i2 = f2iData[ndx2];
    double ut1mutc {f2i1.ut1mutc + dt*(f2i2.ut1mutc - f2i1.ut1mutc)};
    double lod {f2i1.lod + dt*(f2i2.lod - f2i1.lod)};
    Eigen::Quaterniond bpn {f2i1.bpn.slerp(dt, f2i2.bpn)};
    Eigen::Quaterniond pm {f2i1.pm.slerp(dt, f2i2.pm)};
    ecf_eci f2i {mjd2000, ut1mutc, lod, pm, bpn};
    return f2i;
  } else {
    return f2i1;
  }
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::ecf2eci(const JulianDate& utc,
                   const Eigen::Matrix<double, 3, 1>& posf) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(era, Eigen::Vector3d::UnitZ())};
  return f2i.bpn*qera*f2i.pm*posf;
}


Eigen::Matrix<double, 6, 1>
EcfEciSys::ecf2eci(const JulianDate& utc,
                   const Eigen::Matrix<double, 3, 1>& posf,
                   const Eigen::Matrix<double, 3, 1>& velf) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(era, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> pos_tirf = f2i.pm*posf;
  Eigen::Quaterniond qbpn_era {f2i.bpn*qera};
  double we {phy_const::earth_angular_velocity(f2i.lod)};
  Eigen::Matrix<double, 3, 1> wvec {0.0, 0.0, we};

  Eigen::Matrix<double, 3, 1> posi = qbpn_era*pos_tirf;
  Eigen::Matrix<double, 3, 1> veli = qbpn_era*(f2i.pm*velf +
                                               wvec.cross(pos_tirf));

  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3, 1>(0, 0) = posi;
  xeci.block<3, 1>(3, 0) = veli;

  return xeci;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::gravity2ecf(const JulianDate& utc, 
                       const Eigen::Matrix<double, 3, 1>& r_s_o_f,
                       const Eigen::Matrix<double, 3, 1>& v_s_f_f,
                       const Eigen::Matrix<double, 3, 1>& a_s_i_f) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  double we {phy_const::earth_angular_velocity(f2i.lod)};
  Eigen::Matrix<double, 3, 1> wvec {0.0, 0.0, we};

    // Convert position and velocity to instantaneous spin axis
    // coordinates 'w' (Polar motion)
  Eigen::Matrix<double, 3, 1> r_s_o_w = f2i.pm*r_s_o_f;
  Eigen::Matrix<double, 3, 1> v_s_f_w = f2i.pm*v_s_f_f;

  Eigen::Matrix<double, 3, 1> da_s_f_w = 2.0*wvec.cross(v_s_f_w) +
                                             wvec.cross(wvec.cross(r_s_o_w));

  return a_s_i_f - f2i.pm.conjugate()*da_s_f_w;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::eci2ecf(const JulianDate& utc,
                  const Eigen::Matrix<double, 3, 1>& posi) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ())};
  return f2i.pm.conjugate()*qera*f2i.bpn.conjugate()*posi;
}


Eigen::Matrix<double, 6, 1>
EcfEciSys::eci2ecf(const JulianDate& utc,
                   const Eigen::Matrix<double, 3, 1>& posi,
                   const Eigen::Matrix<double, 3, 1>& veli) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ())};
  Eigen::Quaterniond qera_bpnt {qera*f2i.bpn.conjugate()};
  Eigen::Matrix<double, 3, 1> pos_tirf = qera_bpnt*posi;
  double we {phy_const::earth_angular_velocity(f2i.lod)};
  Eigen::Matrix<double, 3, 1> wvec {0.0, 0.0, we};

  Eigen::Matrix<double, 3, 1> posf = f2i.pm.conjugate()*pos_tirf;
  Eigen::Matrix<double, 3, 1> velf = f2i.pm.conjugate()*(qera_bpnt*veli -
                                                         wvec.cross(pos_tirf));
  
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3, 1>(0, 0) = posf;
  xecf.block<3, 1>(3, 0) = velf;

  return xecf;
}


Eigen::Matrix<double, 6, 1>
EcfEciSys::ecf2teme(const JulianDate& utc,
                    const Eigen::Matrix<double, 3, 1>& posf,
                    const Eigen::Matrix<double, 3, 1>& velf) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double gmst {iauGmst82(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qgmst{Eigen::AngleAxisd(gmst, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> pos_tirf = f2i.pm*posf;
  double we {phy_const::earth_angular_velocity(f2i.lod)};
  Eigen::Matrix<double, 3, 1> wvec {0.0, 0.0, we};

  Eigen::Matrix<double, 3, 1> posi = qgmst*pos_tirf;
  Eigen::Matrix<double, 3, 1> veli = qgmst*(f2i.pm*velf + wvec.cross(pos_tirf));

  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3, 1>(0, 0) = posi;
  xeci.block<3, 1>(3, 0) = veli;

  return xeci;
}


Eigen::Matrix<double, 6, 1>
EcfEciSys::teme2ecf(const JulianDate& utc,
                    const Eigen::Matrix<double, 3, 1>& posi,
                    const Eigen::Matrix<double, 3, 1>& veli) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double gmst {iauGmst82(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qgmst{Eigen::AngleAxisd(-gmst, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> pos_tirf = qgmst*posi;
  double we {phy_const::earth_angular_velocity(f2i.lod)};
  Eigen::Matrix<double, 3, 1> wvec {0.0, 0.0, we};

  Eigen::Matrix<double, 3, 1> posf = f2i.pm.conjugate()*pos_tirf;
  Eigen::Matrix<double, 3, 1> velf = f2i.pm.conjugate()*(qgmst*veli -
                                                         wvec.cross(pos_tirf));
  
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3, 1>(0, 0) = posf;
  xecf.block<3, 1>(3, 0) = velf;

  return xecf;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::teme2ecf(const JulianDate& utc,
                    const Eigen::Matrix<double, 3, 1>& posi) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + phy_const::day_per_tu*f2i.ut1mutc};
  double gmst {iauGmst82(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qgmst{Eigen::AngleAxisd(-gmst, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> pos_tirf = qgmst*posi;

  return f2i.pm.conjugate()*pos_tirf;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::mod2eci(const JulianDate& utc,
                   const Eigen::Matrix<double, 3, 1>& mod) const
{
    // Check for valid date
  double days {utc - jdStart};
  if (days < 0.0  ||  jdStop - utc < 0.0) {
    throw std::out_of_range("EcfEciSys::mod2eci() Time out of range");
  }

    // Always need first index
  unsigned long int ndx1 {static_cast<unsigned long int>(days/rate_days)};
  const meme_eci& i2i1 = memeData[ndx1];

    // Get second data set if interpolating - otherwise,
    // return data less than or equal to requested time
  if (nfi == 1UL) {
    return bt*memeData[0UL].p76*mod;
  } else if (interpolate_bpnpm) {
    double mjd2000 {utc.getMjd2000()};
    double dt_days {mjd2000 - i2i1.mjd2000};
    double dt {dt_days/rate_days};
    unsigned long int ndx2 {ndx1 + 1UL};
    const meme_eci& i2i2 = memeData[ndx2];
    Eigen::Quaterniond p76 {i2i1.p76.slerp(dt, i2i2.p76)};
    return bt*p76*mod;
  } else {
    return bt*i2i1.p76*mod;
  }
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::j20002gcrf(const Eigen::Matrix<double, 3, 1>& j2000) const
{
  return bt*j2000;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::gcrf2j2000(const Eigen::Matrix<double, 3, 1>& gcrf) const
{
  return bt.conjugate()*gcrf;
}


}

static Eigen::Matrix3d from3x3(double mtx[3][3])
{
  Eigen::Matrix3d emtx;
  emtx(0,0) = mtx[0][0];
  emtx(1,0) = mtx[1][0];
  emtx(2,0) = mtx[2][0];
  emtx(0,1) = mtx[0][1];
  emtx(1,1) = mtx[1][1];
  emtx(2,1) = mtx[2][1];
  emtx(0,2) = mtx[0][2];
  emtx(1,2) = mtx[1][2];
  emtx(2,2) = mtx[2][2];

  return emtx;
}
