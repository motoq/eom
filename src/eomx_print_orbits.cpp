/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <Eigen/Dense>

#include <astro_ecfeci_sys.h>
#include <astro_keplerian.h>

#include <eomx.h>


/**
 * See eomx.h
 */
void eomx_print_orbits(const std::unordered_map<std::string,
                             std::shared_ptr<eom::Ephemeris>>& ephemerides,
                       const std::string& file_name,
                       const std::shared_ptr<const eom::EcfEciSys>& f2iSys)
{
  std::cout << "\n\nGenerated Orbits";
  for (const auto& [name, eph] : ephemerides) {
    std::cout << "\n  " << name;
    std::cout << "\n  " << eph->getEpoch().to_string() << "    GCRF";
    eom::Keplerian oeCart(eph->getStateVector(eph->getEpoch(),
                                              eom::EphemFrame::eci));
    std::cout << oeCart;
  }

  if (file_name.length() > 0) {
    std::ofstream fout(file_name);
    if (fout.is_open()) {
      fout << "ITRF";
      fout << "\n----";
      for (const auto& [name, eph] : ephemerides) {
        fout << "\n  " << name;
        fout << "\n  " << eph->getEpoch().to_string() << "    ITRF";
        Eigen::Matrix<double, 6, 1> xef =
            eph->getStateVector(eph->getEpoch(), eom::EphemFrame::ecf);
        fout << xef;
      }

      fout << "\nGCRF";
      fout << "\n----";
      for (const auto& [name, eph] : ephemerides) {
        fout << "\n  " << name;
        fout << "\n  " << eph->getEpoch().to_string();
        eom::Keplerian oeCart(eph->getStateVector(eph->getEpoch(),
                                                  eom::EphemFrame::eci));
        fout << oeCart;
      }

      fout << "\nJ2000";
      fout << "\n----";
      for (const auto& [name, eph] : ephemerides) {
        fout << "\n  " << name;
        fout << "\n  " << eph->getEpoch().to_string();
        Eigen::Matrix<double, 6, 1> xeci = 
            eph->getStateVector(eph->getEpoch(), eom::EphemFrame::eci);
        xeci.block<3,1>(0,0) = f2iSys->gcrf2j2000(xeci.block<3,1>(0,0));
        xeci.block<3,1>(3,0) = f2iSys->gcrf2j2000(xeci.block<3,1>(3,0));
        eom::Keplerian oeCart(xeci);
        fout << oeCart;
      }


      fout << '\n';
    } else {
      std::cerr << "\n\n  Invalid Orbit Summary Filename:  " <<
                   file_name << '\n';
    }
  }
}

