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
#include <vector>

#include <Eigen/Dense>

#include <phy_const.h>
#include <astro_composite_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_generate.h>
#include <astro_keplerian.h>

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

  m_xferStart = parse_datetime(tokens);
  m_xferDur = parse_duration(tokens);

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
  auto xferEnd = m_xferStart + m_xferDur;
  if (m_xferStart < m_start_eph->getBeginTime() ||
      m_start_eph->getEndTime() < xferEnd) {
    throw CmdValidateException("EomXferOrbit::validate() "s +
                               m_start_eph->getName() +
                               " does not cover transfer orbit duration"s);
  }
  if (m_xferStart < m_end_eph->getBeginTime() ||
      m_end_eph->getEndTime() < xferEnd) {
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
      m_xferStart,
      m_xferDur,
      m_propCfg,
      m_f2i);

  auto file_name = m_func_name + ".m"s;
  std::ofstream fout(file_name);
  if (fout.is_open()) {
    if (nitr < 0) {
      fout << "\n%  No solution found for transfer orbit\n";
      return;
    }
    fout  << "\n%  Transfer converged in " << nitr << " iterations";

    auto xferEnd = m_xferStart + m_xferDur;
    Eigen::Matrix<double, 6, 1> x1 =
        m_start_eph->getStateVector(m_xferStart, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x1t =
        xfer_eph->getStateVector(m_xferStart, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x2 =
        m_end_eph->getStateVector(xferEnd, eom::EphemFrame::eci);
    Eigen::Matrix<double, 6, 1> x2t =
        xfer_eph->getStateVector(xferEnd, eom::EphemFrame::eci);
    Eigen::Matrix<double, 3, 1> dv1 = x1t.block<3,1>(3,0) -
                                       x1.block<3,1>(3,0);
    Eigen::Matrix<double, 3, 1> dv2 = x2t.block<3,1>(3,0) -
                                       x2.block<3,1>(3,0);

    try {
      fout << ("\n%{");
      eom::Keplerian kep {x1t};
      fout << "\nTransfer Orbit:\n" << kep;
      fout << ("\n%}");
    } catch (const std::invalid_argument& ia) {
      fout << "\n\nCan't form elliptical orbit:  " << ia.what();
    }

    fout << "\n%  Entry DeltaV:  " << phy_const::m_per_du*
                                      phy_const::tu_per_sec*dv1.norm() <<
                                      " m/sec";
    fout << "\n%  Exit DeltaV:   " << phy_const::m_per_du*
                                      phy_const::tu_per_sec*dv2.norm() <<
                                      " m/sec";

      // Composite ephemeris for plotting
    std::vector<eom::JulianDate> ho_times = {m_xferStart, xferEnd};
    std::vector<std::shared_ptr<eom::Ephemeris>> veph = {m_start_eph,
                                                         std::move(xfer_eph),
                                                         m_end_eph};
    eom::CompositeEphemeris ceph {"composite", ho_times, veph};

    auto jdStart = m_xferStart - m_xferDur;
    auto jdStop = xferEnd + m_xferDur;
    double tot_time {m_to_time_units*phy_const::tu_per_day*(jdStop - jdStart)};
    double dt {m_to_time_units*m_dtOut.getTu()};
    long int nrec {static_cast<long int>(tot_time/dt) + 1L};

      // Function header
    fout << "function [gxh, tpv] = " << m_func_name;
    fout << "\n% Orbit is an EOM generated Matlab/Octave function that";
    fout << "\n% plots a 3D orbit trace in ECI coordinates based on";
    fout << "\n% composite ephemerides";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh  Graphics handle to new image";
    fout << "\n%   tpv  Nx7 matrix of time, pos, vel, in ECI";
    fout << "\n%        coordinates, units of DU and DU/TU";
    fout << '\n';

      // Create time and range data
    fout << "\ntpv = [";
    fout << std::scientific;
    fout.precision(16);
    for (long int ii=0L; ii<nrec; ++ii) {
      if (ii > 0L) {
        fout << ';';
      }
      double dtnow {ii*dt};
      fout << "\n  " << dtnow << " ";
      eom::JulianDate jdNow {jdStart +
                             phy_const::day_per_tu*(dtnow/m_to_time_units)};
      std::cout << "   "  << jdNow.to_string();
      Eigen::Matrix<double, 6, 1> pv =
          ceph.getStateVector(jdNow, eom::EphemFrame::eci);
      for (int jj=0; jj<6; ++jj) {
        fout << " " << pv(jj);
      }
    }

      // Make the plot and annotate
    std::string coords = "ECI";
    fout << "\n];";
    fout << "\nn = size(tpv,1);";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot3(tpv(:,2), tpv(:,3), tpv(:,4), '.b');";
    fout << "\nscatter3(tpv(1,2), tpv(1,3), tpv(1,4), 'g');";
    fout << "\nscatter3(tpv(n,2), tpv(n,3), tpv(n,4), 'r');";
    fout << "\nscatter3(0, 0, 0, 'b');";
    fout << "\nplot3(tpv(:,5), tpv(:,6), tpv(:,7), '.m');";
    fout << "\nscatter3(tpv(1,5), tpv(1,6), tpv(1,7), 'g');";
    fout << "\nscatter3(tpv(n,5), tpv(n,6), tpv(n,7), 'r');";
    fout << "\nxlabel('X, dX/dT');";
    fout << "\nylabel('Y, dY/dT');";
    fout << "\nzlabel('Z, dZ/dT');";
    fout << "\ntitle('" << coords << " " << ceph.getName() <<
            " on " << jdStart.to_dmy_str() << "');";
    fout << "\naxis equal;";
    fout << "\nend\n";
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
