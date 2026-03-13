/*
 * Copyright 2026 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_f2i_msg_cmd.h>

#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>

#include <eom_config.h>
#include <eom_parse.h>

/* local utility */
static void write_dcm(std::ofstream& fout, const Eigen::Quaterniond& q);

namespace eom_app {


EomF2iMsgCmd::EomF2iMsgCmd(std::deque<std::string>& tokens,
                           const EomConfig&)
{
  using namespace std::string_literals;
  if (tokens.size() < 2) {
    throw std::invalid_argument("EomF2iMsgCmd::EomF2iMsgCmd() "s +
                                " invalid number of arguments ");
  }

    // let potential invalid_argument pass through
  m_time = parse_datetime(tokens);
  m_filename = tokens[0];
  tokens.pop_front();
}


/**
 * Validation delegated
 */
void EomF2iMsgCmd::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>&)
{
}


/*
 * Get a copy of the ECF2ECI service
 */
void EomF2iMsgCmd::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>&,
    const std::vector<eom::OrbitDef>&,
    std::shared_ptr<const eom::EcfEciSys> ecfeciSys)
{
  using namespace std::string_literals;

  m_f2i = std::move(ecfeciSys);

  if (m_time < m_f2i->getBeginTime()  ||  m_f2i->getEndTime() < m_time) {
    throw CmdValidateException(
        "EomF2iMsgCmd::validate()"s +
        " Requested time not within ECF2ECI data range");
  }
}


/*
 * Open and write ECFECI data to a file for the requested time
 */
void EomF2iMsgCmd::execute() const
{
  using namespace std::string_literals;

  std::ofstream fout(m_filename);
  if (fout.is_open()) {
    fout << '\n' << m_time;
    eom::ecf_eci f2i_msg = m_f2i->getEcfEciData(m_time);
    fout << "\nIAU 2000A CIO ITRF to TIRS (W)";
    write_dcm(fout, f2i_msg.pm);
    fout << "\nIAU 2000A CIO TIRS to CIRS (ERA)";
    fout << "\n    " << utl_const::deg_per_rad*m_f2i->getEra(m_time) << " deg";
    fout << "\nIAU 2000A CIO CIRS to GCRF (BPN)";
    write_dcm(fout, f2i_msg.bpn);
    fout << '\n';
  } else {
    std::cerr << "\nCan't open " << m_filename << '\n';
  }


}

}

static void write_dcm(std::ofstream& fout, const Eigen::Quaterniond& q) {
  Eigen::Matrix<double, 3, 3> dcm(q);
  fout << std::scientific;
  fout.precision(16);
  for (int ii=0; ii<3; ++ii) {
    fout << '\n';
    for (int jj=0; jj<3; ++jj) {
      fout << "    " << dcm(ii,jj);
    }
  }
}
