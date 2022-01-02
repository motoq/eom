#include <astro_keplerian.h>

#include <ostream>
#include <iomanip>
#include <cmath>
#include <array>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <utl_const.h>
#include <phy_const.h>

//namespace {
  //constexpr int    niter {100};
  //constexpr double eps {1.e-10};

  //constexpr int itmax {500};
  //constexpr double ckep {1.0e-20};

  //constexpr double gm {3600.*398600.4418/(phy_const::km_per_du*
  //                                        phy_const::km_per_du*
  //                                        phy_const::km_per_du)};
//}

namespace eom {


/*
 * Based on Vallado's 4th edition, Algorithm 10: COE2RV
 */
Keplerian::Keplerian(const std::array<double, 6>& oe)
{
  m_oe = oe;

  double a {oe[0]};
  double e {oe[1]};
  double i {oe[2]};
  double o {oe[3]};
  double w {oe[4]};
  double v {oe[5]};

  double semip {a*(1.0 - e*e)};
  double cv {std::cos(v)};
  double sv {std::sin(v)};
  double ecv {e*cv};
  double suop {std::sqrt(phy_const::gm/semip)};

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

  m_xyz.block<3,1>(0,0) = q_pqw2eci*r_pqw;
  m_xyz.block<3,1>(3,0) = q_pqw2eci*v_pqw;
}


/*
 * Based on Vallado's 4th edition, Algorithm 9: RV2COE
 */
Keplerian::Keplerian(const Eigen::Matrix<double, 6, 1>& xyz)
{
  m_xyz = xyz;

  Eigen::Matrix<double, 3, 1> rvec {m_xyz.block<3,1>(0,0)};
  Eigen::Matrix<double, 3, 1> vvec {m_xyz.block<3,1>(3,0)};

  Eigen::Matrix<double, 3, 1> hvec {rvec.cross(vvec)};
  Eigen::Matrix<double, 3, 1> khat {Eigen::Vector3d::UnitZ()};
  Eigen::Matrix<double, 3, 1> nvec {khat.cross(hvec)};

  double rmag {rvec.norm()};
  double vmag {vvec.norm()};
  double hmag {hvec.norm()};
  double nmag {nvec.norm()};
  double v2 {vmag*vmag};
  double rdotv {rvec.dot(vvec)};
  double muor {phy_const::gm/rmag};

    // Eccentricity
  Eigen::Matrix<double, 3, 1> evec {((v2 - muor)*rvec - rdotv*vvec)/
                                     phy_const::gm
                                   };
  double emag {evec.norm()};
    // Vis-viva eqn
  double sme {v2/2.0 - muor};
    // Semimajor axis
  double sma {-0.5*phy_const::gm/sme};
    // Inclination
  double inc {std::acos(hvec(2)/hmag)};
    // RAAN
  double raan {std::acos(nvec(0)/nmag)};
  if (nvec(1) < 0.0) {
    raan = utl_const::tpi - raan;
  }
    // Argument of perigee
  Eigen::Matrix<double, 3, 1> ehat {evec/emag};
  Eigen::Matrix<double, 3, 1> nhat {nvec/nmag};
  double argp {std::acos(nhat.dot(ehat))};
  if (evec(2) < 0.0) {
    argp = utl_const::tpi - argp;
  }
    // True anomaly
  Eigen::Matrix<double, 3, 1> rhat {rvec/rmag};
  double ta {std::acos(ehat.dot(rhat))};
  if (rdotv < 0.0) {
    ta = utl_const::tpi - ta;
  }

  m_oe[0] = sma;
  m_oe[1] = emag;
  m_oe[2] = inc;
  m_oe[3] = raan;
  m_oe[4] = argp;
  m_oe[5] = ta;
}


/*
double trueAnomaly(double ma, double e)
{
  return ta;
}


double meanAnomaly(double f, double e)
{
  return ea - e*sea;
}
*/

void Keplerian::print(std::ostream& stream) const
{
  stream << std::fixed;
  stream << std::setprecision(3);
  stream << "\n  a:" << phy_const::km_per_du*m_oe[0] << " km";
  stream << std::setprecision(5);
  stream << " e:" << m_oe[1];
  stream << " i:" << utl_const::deg_per_rad*m_oe[2] << "\u00B0";
  stream << " o:" << utl_const::deg_per_rad*m_oe[3] << "\u00B0";
  stream << " w:" << utl_const::deg_per_rad*m_oe[4] << "\u00B0";
  stream << " v:" << utl_const::deg_per_rad*m_oe[5] << "\u00B0";
  stream << std::setprecision(3);
  stream << "\n    {" << phy_const::km_per_du*m_xyz(0) << "  " <<
                         phy_const::km_per_du*m_xyz(1) << "  " <<
                         phy_const::km_per_du*m_xyz(2) << "} km";
  stream << std::setprecision(6);
  stream << "\n    {" <<
            phy_const::km_per_du*m_xyz(3)*phy_const::tu_per_sec << "  " <<
            phy_const::km_per_du*m_xyz(4)*phy_const::tu_per_sec << "  " <<
            phy_const::km_per_du*m_xyz(5)*phy_const::tu_per_sec << "} km/sec";
}


}

