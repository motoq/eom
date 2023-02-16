/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <string>
#include <memory>

#include <Eigen/Dense>

#include <cal_duration.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_moon_meeus.h>
#include <astro_print.h>

#include <eom_test.h>

namespace eom_app {

void eom_test_moon()
{
  std::cout << "\n\n  === Test:  Moon ===";
  std::cout << "\n  Generating Meeus based moon ephemeris";

  /*
   * Meeus Example:eom::GregDate gdStart(1992, 10, 13);
  eom::GregDate gdStart(1992, 4, 12);
  double tott {-1.0*(37.0 + 32.184)/86400.0};
  eom::JulianDate jdStart(gdStart);
  jdStart += tott;
   */

  eom::GregDate gdStart(2023, 02, 04);
  eom::JulianDate jdStart(gdStart);
  auto jdStop = jdStart + 30.0;

  eom::Duration dt {1.0, phy_const::tu_per_day};
  std::shared_ptr<eom::EcfEciSys> ecfeci =
          std::make_shared<eom::EcfEciSys>(jdStart,
                                           jdStop+1.0,
                                           dt,
                                           nullptr);

  std::shared_ptr<const eom::Ephemeris> ephPtr =
          std::make_shared<eom::MoonMeeus>(ecfeci);


  eom::Duration dtEph {1.25, phy_const::tu_per_day};
  eom::print_ephemeris("MoonMeeus.e", jdStart, jdStop, dtEph,
                       eom::EphemFrame::eci, ephPtr);


  std::cout << "\n  === End Test:  Moon ===\n\n";
}


}

