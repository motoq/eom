/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_gravity_std.h>

#include <iostream>
#include <cmath>

#include <string>
#include <memory>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <mth_ode.h>
#include <astro_egm_coeff.h>
#include <astro_gravity.h>
#include <mth_legendre_af.h>

namespace eom {

/*
 * Coefficients in astro_egm_coeff.h are expected to be unormalized and
 * continuous, sorted by order first, then degree.
 */
GravityStd::GravityStd(int max_degree, int max_order)
{
  if (max_order > max_degree) {
    throw std::invalid_argument("GravityStd::GravityStd Order > Degree");
  }
  if (max_degree < 0  ||  max_degree > egm_coeff::degree) {
    throw std::invalid_argument(
        "GravityStd::GravityStd() Unsupported Degree: " +
        std::to_string(max_degree)
    );
  }
  if (max_order < 0  ||  max_order > egm_coeff::order) {
    throw std::invalid_argument(
        "GravityStd::GravityStd() Unsupported Order: " +
        std::to_string(max_order)
    );
  }

  m_degree = max_degree;
  m_order = max_order;
  m_smlon = std::make_unique<double[]>(max_order + 1);
  m_cmlon = std::make_unique<double[]>(max_order + 1);
  m_re_r_n = std::make_unique<double[]>(max_degree + 1);
  m_alf = std::make_unique<LegendreAf>(m_degree, m_order);
}


Eigen::Matrix<double, 3, 1>
    GravityStd::getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                           OdeEvalMethod entry)
{
  using namespace egm_coeff;

    // Geometry
  const double rx {pos(0)};
  const double ry {pos(1)};
  const double rz {pos(2)};
  const double rmag {pos.norm()};
  const double invr {1.0/rmag};
  const double invr2 {invr*invr};
  const double rxy2 {rx*rx + ry*ry};
  const double rxy {std::sqrt(rxy2)};
  const double invrxy {1.0/rxy};
  const double invrxy2 {invrxy*invrxy};
  const double slat {rz*invr};
  const double clat {rxy*invr};
  const double tlat {rz*invrxy};
  const double slon {ry*invrxy};
  const double clon {rx*invrxy};
  const double re_r {phy_const::re*invr};
  const double gm_r {phy_const::gm*invr};

  m_re_r_n[0] = 1.0;
  m_re_r_n[1] = re_r;
  m_smlon[0] = 0.0;
  m_smlon[1] = slon;
  m_cmlon[0] = 1.0;
  m_cmlon[1] = clon;
  for (int mdx=2; mdx<=m_order; ++mdx) {
    m_smlon[mdx] = 2*clon*m_smlon[mdx-1] - m_smlon[mdx-2];
    m_cmlon[mdx] = 2*clon*m_cmlon[mdx-1] - m_cmlon[mdx-2];
    m_re_r_n[mdx] = re_r*m_re_r_n[mdx-1];
  }
  for (int ndx=(m_order+1); ndx<=m_degree; ++ndx) {
    m_re_r_n[ndx] = re_r*m_re_r_n[ndx-1];
  }

  m_alf->set(slat, clat);

    // Init accumulation of partials w.r.t. spherical
  double du_dr {0.0};
  double du_dlat {0.0};
  double du_dlon {0.0};
    // Loop based on indexing in egm_coeff
  for (int ndx=0; ndx<nc; ++ndx) {
    const auto nn = xn[ndx];
    const auto mm = xm[ndx];
    if (nn <= m_degree  &&  mm <= m_order) {
      double smlon {m_smlon[mm]};
      double cmlon {m_cmlon[mm]};
      double pnm = (*m_alf)(nn, mm);
      double pnmp1 = (*m_alf)(nn, mm+1);
      {
        double lon = atan2(slon, clon);
        if (std::fabs(smlon - std::sin(mm*lon)) > 1.0e-12  ||
            std::fabs(cmlon - std::cos(mm*lon)) > 1.0e-12) {
          std::cerr << "\nRecursive trig error\n";
        }
      }
      double scnm = cnm[ndx];
      double ssnm = snm[ndx];
      du_dr += (nn+1)*m_re_r_n[nn]*pnm*(scnm*cmlon + ssnm*smlon);
      du_dlat += m_re_r_n[nn]*(pnmp1 - mm*tlat*pnm)*
                              (scnm*cmlon + ssnm*smlon);
      du_dlon += mm*m_re_r_n[nn]*pnm*(ssnm*cmlon - scnm*smlon);
/*
      du_dr += (nn + 1)*m_re_r_n[nn]*(*m_alf)(nn, mm)*
                       (cnm[ndx]*m_cmlon[mm] + snm[ndx]*m_smlon[mm]);
      du_dlat += m_re_r_n[nn]*((*m_alf)(nn, mm+1) - mm*tlat*(*m_alf)(nn, mm))*
                              (cnm[ndx]*m_cmlon[mm] + snm[ndx]*m_smlon[mm]);
      du_dlon += mm*m_re_r_n[nn]*(*m_alf)(nn, mm)*
                                 (snm[ndx]*m_cmlon[mm] - cnm[ndx]*m_smlon[mm]);
*/
    }
    ndx++;
  }
    // Finish partials
  du_dr += 1.0;
  du_dr *= -1.0*gm_r*invr;
  du_dlat *= gm_r;
  du_dlon *= gm_r;
    // Convert from spherical to Cartesian
  double dlat {invr*du_dr - du_dlat*rz*invrxy*invr2};
  double dlon {du_dlon*invrxy2};
  double ax {dlat*rx - dlon*ry};
  double ay {dlat*ry + dlon*rx};
  double az {invr*du_dr*rz + du_dlat*rxy*invr2};

  Eigen::Matrix<double, 3, 1> acc = {ax, ay, az};
  
  return acc;
  //return -1.0*phy_const::gm*invr*invr*invr*pos;
}


}
