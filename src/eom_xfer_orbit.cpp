/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_xfer_orbit.h>

#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <Eigen/Dense>
/*


#include <cal_julian_date.h>
#include <astro_attitude.h>
*/

#include <phy_const.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_generate.h>

#include <eom_config.h>
#include <eom_parse.h>

namespace eom_app {


EomXferOrbit::EomXferOrbit(std::deque<std::string>& tokens,
                           const EomConfig& cfg)
{
    // Read orbit name, etc.
  using namespace std::string_literals;
  if (tokens.size() != 12) {
    throw std::invalid_argument("EomXferOrbit::EomXferOrbit() "s +
                                "TransferOrbit requires 12 arguments "s +
                                "vs. input "s +
                                std::to_string(tokens.size()));
  }
  m_start_orbit_name = tokens[0];
  tokens.pop_front();
  if (!cfg.pendingOrbit(m_start_orbit_name)) {
      throw std::invalid_argument("EomXferOrbit::EomXferOrbit() "s +
                                  "Invalid start orbit " + m_start_orbit_name);
  }
  m_end_orbit_name = tokens[0];
  tokens.pop_front();
  if (!cfg.pendingOrbit(m_end_orbit_name)) {
      throw std::invalid_argument("EomXferOrbit::EomXferOrbit() "s +
                                  "Invalid end orbit " + m_end_orbit_name);
  }

  m_xfer_start = parse_datetime(tokens);
  m_xfer_dur = parse_duration(tokens);

  m_func_name = tokens[0];
  tokens.pop_front();

  m_dtOut = cfg.getOutputRate();
  m_distanceUnitsLbl = cfg.getIoDistansUnits();
  m_to_time_units = cfg.getIoPerTu();
  m_to_distance_units = cfg.getIoPerDu();
}


/**
 * Validation delegated
 */
void EomXferOrbit::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>&)
{
}


/*
 * Set ephemeris pointers using orbit_names from initialization
 */
void EomXferOrbit::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides,
    const std::vector<eom::OrbitDef>& orbits,
    std::shared_ptr<const eom::EcfEciSys> ecfeciSys)
{
  using namespace std::string_literals;

  m_f2i = std::move(ecfeciSys);

    // Connect with ephemeris resources
  try {
    m_start_eph = ephemerides.at(m_start_orbit_name);
    m_end_eph = ephemerides.at(m_end_orbit_name);
  } catch (const std::out_of_range& oor) {
    throw CmdValidateException("EomXferOrbit::validate() "s +
                               "Invalid orbit name in TransferOrbit Cmd: "s +
                               m_start_orbit_name + " or "s + m_end_orbit_name);
  }

    // Use orbit propagator for start orbit to reach end orbit state
  bool found_propagator {false};
  for (const auto& odef : orbits) {
    if (odef.getOrbitName() == m_start_eph->getName()) {
      m_propCfg = odef.getPropagatorConfig();
      found_propagator = true;
      break;
    }
  }
  if (!found_propagator) {
    throw CmdValidateException(
        "EomXferOrbit::validate()"s +
        "\n  Propagator not found in TransferOrbit Cmd for: "s +
        m_start_eph->getName() +
        "\n  Ensure start orbit is neither a derived orbit nor ephemeris");
  }

    // Check that supplied ephemerides cover transfer orbit time span
  auto xfer_end = m_xfer_start + m_xfer_dur;
  if (m_xfer_start < m_start_eph->getBeginTime() ||
      m_start_eph->getEndTime() < xfer_end) {
    throw CmdValidateException("EomXferOrbit::validate() "s +
                               m_start_eph->getName() +
                               " does not cover transfer orbit duration"s);
  }
  if (m_xfer_start < m_end_eph->getBeginTime() ||
      m_end_eph->getEndTime() < xfer_end) {
    throw CmdValidateException("EomXferOrbit::validate() "s +
                               m_end_eph->getName() +
                               " does not cover transfer orbit duration"s);
  }

}


void EomXferOrbit::execute() const
{
  using namespace std::string_literals;

  //std::unique_ptr<eom::Ephemeris> xfer_eph = eom::generate_xfer_orbit(
  auto [xfer_eph, nitr] = eom::generate_xfer_orbit(
      m_start_orbit_name + "_to_"s + m_end_orbit_name,
      *m_start_eph,
      *m_end_eph,
      m_xfer_start,
      m_xfer_dur,
      m_propCfg,
      m_f2i);

  auto file_name = m_func_name + ".lis"s;
  std::ofstream fout(file_name);
  if (fout.is_open()) {
    if (nitr < 0) {
      fout << "\nNo solution found for transfer orbit\n";
      return;
    }
    auto xfer_end = m_xfer_start + m_xfer_dur;

    fout  << "\nTransfer converged in " << nitr << " iterations";
    Eigen::Matrix<double, 6, 1> x1 =
        m_start_eph->getStateVector(m_xfer_start, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x1t =
        xfer_eph->getStateVector(m_xfer_start, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x2 =
        m_end_eph->getStateVector(xfer_end, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x2t =
        xfer_eph->getStateVector(xfer_end, eom::EphemFrame::eci);
    Eigen::Matrix<double, 3, 1> dv1 = x1t.block<3,1>(3,0) -
                                       x1.block<3,1>(3,0);
    Eigen::Matrix<double, 3, 1> dv2 = x2t.block<3,1>(3,0) -
                                       x2.block<3,1>(3,0);
    fout << "\nEntry DeltaV:  " << phy_const::m_per_du*
                                   phy_const::tu_per_sec*dv1.norm() <<
                                        " m/sec";
    fout << "\nExit DeltaV:   " << phy_const::m_per_du*
                                   phy_const::tu_per_sec*dv2.norm() << " m/sec";

    



/*
    double tot_time {m_to_time_units*phy_const::tu_per_day*(m_jdStop -
                                                            m_jdStart)};
    double dt {m_to_time_units*m_dtOut.getTu()};
    unsigned long int nrec {static_cast<unsigned long int>(tot_time/dt)};
    nrec++;

      // Function header
    fout << "function [gxh, time_rtc] = " << m_func_name;
    fout << "\n% RTC is an EOM generated Matlab/Octave function that";
    fout << "\n% plots relative position as a function of time";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh       Graphics handle to new image";
    fout << "\n%   time_rtc  Nx4 matrix of time and RTC coordinates";
    fout << '\n';

      // Create time and range data
    fout << "\ntime_rtc = [";
    fout << std::scientific;
    fout.precision(16);
    for (unsigned long int ii=0UL; ii<nrec; ++ii) {
      if (ii > 0UL) {
        fout << ';';
      }
      double dtnow {ii*dt};
      fout << "\n  " << dtnow << " ";
      eom::JulianDate jdNow {m_jdStart +
                             phy_const::day_per_tu*(dtnow/m_to_time_units)};
      Eigen::Matrix<double, 6, 1> pv1 =
          m_eph[0]->getStateVector(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> r1 {pv1.block<3,1>(0,0)};
      Eigen::Matrix<double, 3, 1> r2 =
          m_eph[1]->getPosition(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> dr = r1 - r2;
      Eigen::Matrix<double, 3, 1> v1 {pv1.block<3,1>(3,0)};
      Eigen::Matrix<double, 3, 3> i2rtcDcm {eom::AttitudeRtc<double>(r1, v1)};
      dr = i2rtcDcm*dr;
      for (int jj=0; jj<3; ++jj) {
        fout << " " << m_to_distance_units*dr(jj);
      }
    }
      // Make the plot and annotate
    fout << "\n];";
    fout << "\nn = size(time_rtc,1);";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot3(time_rtc(:,2), time_rtc(:,3), time_rtc(:,4));";
    fout << "\nscatter3(time_rtc(1,2), time_rtc(1,3), time_rtc(1,4), 'g');";
    fout << "\nscatter3(time_rtc(n,2), time_rtc(n,3), time_rtc(n,4), 'r');";
    fout << "\nscatter3(0, 0, 0, 'b');";
    fout << "\nxlabel('Radial (" << m_distanceUnitsLbl << ")');";
    fout << "\nylabel('Transverse (" << m_distanceUnitsLbl << ")');";
    fout << "\nzlabel('Cross-Track (" << m_distanceUnitsLbl << ")');";
    fout << "\ntitle('" << m_orbit_names[0] << '-' << m_orbit_names[1] <<
            " RTC on " << m_jdStart.to_dmy_str() << "');";
    fout << "\naxis equal;";
    fout << "\nend\n";
    fout.close();
*/
    fout << '\n';
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
