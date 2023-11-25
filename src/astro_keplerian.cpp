/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_keplerian.h>

#include <ostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <utl_nonconvergence_exception.h>
#include <utl_const.h>
#include <phy_const.h>
#include <mth_angle.h>

namespace {
    // Gravitational parameter
  constexpr double gm {phy_const::gm};
    // Indexing
  constexpr int ia {0};           // Semimajor axis
  constexpr int ie {1};           // Eccentricity
  constexpr int ii {2};           // Inclination
  constexpr int io {3};           // RAAN
  constexpr int iw {4};           // Argument of perigee
  constexpr int iv {5};           // True anomaly
    // Convergence
  constexpr int niter {100};
  constexpr double eps {1.e-10};
    // oe_eps {
  constexpr double oe_eps {1.0e-6};
}

namespace eom {

Keplerian::Keplerian(const std::array<double, 6>& oe)
{
  this->set(oe);
}

/*
 * Based on Vallado's "Fundamentals of Astrodynamics and Applications",
 * 4th edition, Algorithm 9: RV2COE
 */
Keplerian::Keplerian(const Eigen::Matrix<double, 6, 1>& cart)
{
  m_cart = cart;

  Eigen::Matrix<double, 3, 1> rvec {m_cart.block<3,1>(0,0)};
  Eigen::Matrix<double, 3, 1> vvec {m_cart.block<3,1>(3,0)};

  Eigen::Matrix<double, 3, 1> hvec {rvec.cross(vvec)};
  Eigen::Matrix<double, 3, 1> khat {Eigen::Vector3d::UnitZ()};
  Eigen::Matrix<double, 3, 1> nvec {khat.cross(hvec)};
  double nmag {nvec.norm()};
    // By definition - potential numerical roundoff with cross product
  nvec(2) = 0.0;

  double rmag {rvec.norm()};
  double vmag {vvec.norm()};
  m_hmag = hvec.norm();
  double v2 {vmag*vmag};
  double rdotv {rvec.dot(vvec)};
  double muor {gm/rmag};

    // Eccentricity
  Eigen::Matrix<double, 3, 1> evec {((v2 - muor)*rvec - rdotv*vvec)/gm};
  double emag {evec.norm()};
    // Vis-viva eqn
  m_sme = v2/2.0 - muor;

    // Some initial error checking
  if (emag < oe_eps) {
    throw std::invalid_argument(
        "Keplerian::Keplerian(): Eccentricity too close to zero");
  }
  if (nmag < oe_eps) {
    throw std::invalid_argument(
        "Keplerian::Keplerian(): Inclination too close to zero");
  }
  if (m_sme >= 0.0) {
    throw std::invalid_argument(
        "Keplerian::Keplerian(): Orbit must be elliptical");
  }
    // Semimajor axis, perigee radius, final error check
  double sma {-0.5*gm/m_sme};
  m_rp = sma*(1.0 - emag);
  if (m_rp < phy_const::re) {
    throw std::invalid_argument(
        "Keplerian::set(): Perigee distance less than 1 DU");
  }

    // Inclination
  double inc {};
  if (hvec(0) == 0.0  &&  hvec(1) == 0.0) {
    if (hvec(2) > 0.0) {
      inc = 0.0;
    } else {
      inc = utl_const::pi; 
    }
  } else {
    inc = std::acos(hvec(2)/m_hmag);
  }
    // RAAN
  double raan {};
  if (nvec(1) == 0.0) {
    raan = 0.0;
  } else {
    raan = std::acos(nvec(0)/nmag);
  }
  if (nvec(1) < 0.0) {
    raan = utl_const::tpi - raan;
  }
    // Argument of perigee
  Eigen::Matrix<double, 3, 1> ehat {evec/emag};
  Eigen::Matrix<double, 3, 1> nhat {nvec/nmag};
  double argp {mth_angle::unit_vec_angle<double>(nhat, ehat)};
  if (evec(2) < 0.0) {
    argp = utl_const::tpi - argp;
  }
    // True anomaly
  Eigen::Matrix<double, 3, 1> rhat {rvec/rmag};
  double ta {mth_angle::unit_vec_angle<double>(ehat, rhat)};
  if (rdotv < 0.0) {
    ta = utl_const::tpi - ta;
  }
  if (ta >= utl_const::tpi) {
    ta -= utl_const::tpi;
  }

  m_oe[ia] = sma;
  m_oe[ie] = emag;
  m_oe[ii] = inc;
  m_oe[io] = raan;
  m_oe[iw] = argp;
  m_oe[iv] = ta;

}


/*
 * Based on Vallado's "Fundamentals of Astrodynamics and Applications",
 * 4th edition, Algorithm 10: COE2RV
 */
void Keplerian::set(const std::array<double, 6>& oe)
{
  m_oe = oe;

  double a {oe[0]};
  double e {oe[1]};
  double i {oe[2]};
  double o {oe[3]};
  double w {oe[4]};
  double v {oe[5]};
  m_rp = a*(1.0 - e);

    // Error checking.  Gravity models not valid below scaling radius.
    // Etc...
  if (m_rp < phy_const::re) {
    throw std::invalid_argument(
        "Keplerian::set(): Perigee distance less than 1 DU");
  }
  if (e < oe_eps) {
    throw std::invalid_argument(
        "Keplerian::set(): Eccentricity too close to zero or negative");
  }
  if (i < oe_eps) {
    throw std::invalid_argument(
        "Keplerian::set(): Inclination too close to zero");
  }

  double semip {a*(1.0 - e*e)};
  double cv {std::cos(v)};
  double sv {std::sin(v)};
  double ecv {e*cv};
  double suop {std::sqrt(gm/semip)};

  m_sme = -0.5*gm/a;
  m_hmag = std::sqrt(gm*semip);

  Eigen::Matrix<double, 3, 1> r_pqw;
  r_pqw(0,0) = semip*cv/(1.0 + ecv);
  r_pqw(1,0) = semip*sv/(1.0 + ecv);
  r_pqw(2,0) = 0.0;
  Eigen::Matrix<double, 3, 1> v_pqw;
  v_pqw(0,0) = -suop*sv;
  v_pqw(1,0) = suop*(e + cv);
  v_pqw(2,0) = 0.0;

  Eigen::Quaterniond qw{Eigen::AngleAxisd(w, Eigen::Vector3d::UnitZ())};
  Eigen::Quaterniond qi{Eigen::AngleAxisd(i, Eigen::Vector3d::UnitX())};
  Eigen::Quaterniond qo{Eigen::AngleAxisd(o, Eigen::Vector3d::UnitZ())};

  Eigen::Quaterniond q_pqw2eci {qo*qi*qw};

  m_cart.block<3,1>(0,0) = q_pqw2eci*r_pqw;
  m_cart.block<3,1>(3,0) = q_pqw2eci*v_pqw;
}


double Keplerian::getEccentricity() const
{
  return m_oe[ie];
}


/*
 * Based on Vallado's "Fundamentals of Astrodynamics and Applications",
 * 4th edition, Algorithm 5: RV2COE
 */
double Keplerian::getEccentricAnomaly() const
{
  double e {m_oe[ie]};
  double v {m_oe[iv]};
  double denom {1.0/(1.0 + e*std::cos(v))};
  double se {std::sin(v)*std::sqrt(1.0 - e*e)*denom};
  double ce {(e + std::cos(v))*denom};

  return std::atan2(se, ce);
}


double Keplerian::getMeanAnomaly() const
{
  double ea {this->getEccentricAnomaly()};
  return  ea - m_oe[ie]*std::sin(ea);
}


double Keplerian::getPeriod() const
{
  double a {m_oe[ia]};
  return utl_const::tpi*std::sqrt(a*a*a/gm);
}


double Keplerian::getPerigeeSpeed() const
{
  return std::sqrt(gm*(2.0/m_rp - 1.0/m_oe[ia]));
}

/*
 * Based on Vallado's "Fundamentals of Astrodynamics and Applications",
 * 4th edition, Algorithms 2 KepEqtnE and 6 Anomaly to
 */
void Keplerian::setWithMeanAnomaly(double ma)
{
  std::array<double, 6> oe = m_oe;

  double e {oe[ie]};
    // Initial guess to E set to M modified by e
  double ea {ma + e};
  if (ma > utl_const::pi  ||
      ma > -utl_const::pi && ma < 0.0) {
    ea = ma - e;
  }

  int itr {0};
  double ea0 {ea};
  while (itr < niter) {
    ea = ea0 + (ma - ea0 + e*std::sin(ea0))/(1.0 - e*std::cos(ea0));
    if (std::fabs(ea - ea0) < eps) {
      break;
    }
    ea0 = ea;
    itr++;
  }
  if (itr == niter) {
    throw NonconvergenceException("Keplerian::setWithMeanAnomaly()");
  }

  double cea {std::cos(ea)};
  double denom {1.0/(1.0 - e*cea)};
  double sv {std::sin(ea)*std::sqrt(1.0 - e*e)*denom};
  double cv {(cea - e)*denom};

  oe[iv] = std::atan2(sv, cv);
  this->set(oe);
}


void Keplerian::print(std::ostream& stream) const
{
  stream << std::fixed;
  stream << std::setprecision(2);
  stream << "    (" << phy_const::tu_per_day/getPeriod() << " rev/day)";
  stream << std::setprecision(3);
  stream << "\n  a: " << phy_const::km_per_du*m_oe[0] << " km";
  stream << std::setprecision(6);
  stream << "  e: " << m_oe[1];
  stream << "  i: " << utl_const::deg_per_rad*m_oe[2] << "\u00B0";
  stream << "  o: " << utl_const::deg_per_rad*m_oe[3] << "\u00B0";
  stream << "  w: " << utl_const::deg_per_rad*m_oe[4] << "\u00B0";
  stream << "\n  v: " << utl_const::deg_per_rad*m_oe[5] << "\u00B0";
  stream << "  M: " <<
            utl_const::deg_per_rad*getMeanAnomaly() << "\u00B0";
  stream << "  E: " <<
            utl_const::deg_per_rad*getEccentricAnomaly() << "\u00B0";
  stream << std::setprecision(3);
  stream << "\n    {" << phy_const::km_per_du*m_cart(0) << "  " <<
                         phy_const::km_per_du*m_cart(1) << "  " <<
                         phy_const::km_per_du*m_cart(2) << "} km";
  stream << std::setprecision(6);
  stream << "\n    {" <<
            phy_const::km_per_du*m_cart(3)*phy_const::tu_per_sec << "  " <<
            phy_const::km_per_du*m_cart(4)*phy_const::tu_per_sec << "  " <<
            phy_const::km_per_du*m_cart(5)*phy_const::tu_per_sec << "} km/sec";
}


}

