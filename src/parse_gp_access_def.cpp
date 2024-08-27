/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <array>
#include <string>
#include <deque>
#include <utility>
#include <stdexcept>

#include <astro_ground_point.h>
#include <axs_gp_constraints.h>
#include <axs_gp_access_def.h>

#include <eom_config.h>

static void parse_constraints(std::deque<std::string>& cnst_toks,
                              const eom_app::EomConfig& cfg, 
                              eom::GpConstraints& constraints,
                              eom::aux_gp_constraints& aux_constraints);

namespace eom_app {

eom::GpAccessDef
parse_gp_access_def(std::deque<std::string>& tokens, const EomConfig& cfg)
{

  using namespace std::string_literals;
    // Need at least the orbit name, ground point name, and access model
  if (tokens.size() < 3) {
    throw std::invalid_argument("eom_app::parse_gp_access_def() "s +
                                "Minimum  of 3 tokens required vs. " +
                                std::to_string(tokens.size()));
  }
  auto method_name = tokens[0];
  tokens.pop_front();
  eom::AccessModel method = eom::AccessModel::std;
  if (method_name == "Standard") {
    method = eom::AccessModel::std;
  } else if (method_name == "Debug") {
    method = eom::AccessModel::dbg;
  } else {
    throw std::invalid_argument("eom_app::parse_gp_access_def() "s +
                                "Invalid Access Algorithm:  " + method_name);
  }
  auto orbit_name = tokens[0];
  tokens.pop_front();
  auto gp_name = tokens[0];
  tokens.pop_front();

    // Loop over all potential options - only one entry per
    // option allowed otherwise exit before using all tokens
    //
    // 1: Minimum elevation
    // 2: Maximum elevation
    // 3: Minimum and maximum azimuth
    // 4: Maximum sun elevation
  int n_constraints {4};
  eom::GpConstraints xcs;
  eom::aux_gp_constraints axcs;
  for (int ii=0; ii<n_constraints; ++ii) {
    parse_constraints(tokens, cfg, xcs, axcs);
    if (tokens.size() == 0) {
        break;
    }
  }

  eom::GpAccessDef gpDef(orbit_name, gp_name, xcs, axcs, method);
  return gpDef;
}


}


// Add return of constraints structure
// Make sure to update number of iterations for caller to
// parse all options listed below
static void parse_constraints(std::deque<std::string>& cnst_toks,
                              const eom_app::EomConfig& cfg, 
                              eom::GpConstraints& constraints,
                              eom::aux_gp_constraints& aux_constraints)
{
  using namespace std::string_literals;

  double rad_per_io {1.0/cfg.getIoPerRad()};
  if (cnst_toks.size() > 1  &&  cnst_toks[0] == "MinimumElevation") {
    cnst_toks.pop_front();
    try {
      constraints.setMinEl(rad_per_io*std::stod(cnst_toks[0]));
      cnst_toks.pop_front();
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_access_def() "s +
                                  "invalid Minimum Elevation: "s +
                                  ia.what());
    }
  } else if (cnst_toks.size() > 1  &&  cnst_toks[0] == "MaximumElevation") {
    cnst_toks.pop_front();
    try {
      constraints.setMaxEl(rad_per_io*std::stod(cnst_toks[0]));
      cnst_toks.pop_front();
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_access_def() "s +
                                  "invalid Maximum Elevation: "s +
                                  ia.what());
    }
  } else if (cnst_toks.size() > 2  &&
             cnst_toks[0] == "MinimumMaximumAzimuth") {
    cnst_toks.pop_front();
    try {
      constraints.setMinMaxAz(rad_per_io*std::stod(cnst_toks[0]),
                              rad_per_io*std::stod(cnst_toks[1]));
      cnst_toks.pop_front();
      cnst_toks.pop_front();
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_access_def() "s +
                                  "invalid Min or Max Azimuth: "s +
                                  ia.what());
    }
  } else if (cnst_toks.size() > 2  &&
             cnst_toks[0] == "SunConstraint") {
    cnst_toks.pop_front();
    if (cnst_toks[0] == "MaximumElevation") {
      cnst_toks.pop_front();
      try {
        aux_constraints.max_sun_el = rad_per_io*std::stod(cnst_toks[0]);
        aux_constraints.use_max_sun_el = true;
        cnst_toks.pop_front();
      } catch (const std::invalid_argument& ia) {
        throw std::invalid_argument("eom_app::parse_access_def() "s +
                                    "invalid Max Sun El: "s +
                                    ia.what());
      }
    }
  }
}

