#include <astro_keplerian.h>

#include <stdexcept>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <array>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <utl_nonconvergence_exception.h>
#include <utl_const.h>
#include <phy_const.h>

namespace {
    // Indexing
  constexpr int ia {0};
  constexpr int ie {1};
  constexpr int ii {2};
  constexpr int io {3};
  constexpr int iw {4};
  constexpr int iv {5};
    // Convergence
  constexpr int niter {100};
  constexpr double eps {1.e-10};
}

namespace eom {

/*
 * Based on Vallado's 4th edition, Algorithm 10: COE2RV
 */
Keplerian::Keplerian(const std::array<double, 6>& oe)
{
  this->set(oe);
}


void Keplerian::set(const std::array<double, 6>& oe)
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

  hmag = std::sqrt(phy_const::gm*semip);

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


/*
 * Based on Vallado's 4th edition, Algorithm 9: RV2COE
 */
Keplerian::Keplerian(const Eigen::Matrix<double, 6, 1>& cart)
{
  m_cart = cart;

  Eigen::Matrix<double, 3, 1> rvec {m_cart.block<3,1>(0,0)};
  Eigen::Matrix<double, 3, 1> vvec {m_cart.block<3,1>(3,0)};

  Eigen::Matrix<double, 3, 1> hvec {rvec.cross(vvec)};
  Eigen::Matrix<double, 3, 1> khat {Eigen::Vector3d::UnitZ()};
  Eigen::Matrix<double, 3, 1> nvec {khat.cross(hvec)};

  double rmag {rvec.norm()};
  double vmag {vvec.norm()};
  hmag = hvec.norm();
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

  m_oe[ia] = sma;
  m_oe[ie] = emag;
  m_oe[ii] = inc;
  m_oe[io] = raan;
  m_oe[iw] = argp;
  m_oe[iv] = ta;
}


/*
 * Based on Vallado's 4th edition, Algorithm 5: RV2COE
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


/*
 * Based on Vallado's 4th edition, Algorithms 2 KepEqtnE and 6 Anomaly to
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
    throw NonconvergenceException("Keplerian::setWithMeanAnomaly");
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

