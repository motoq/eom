/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <Eigen/Dense>

#include <phy_const.h>
#include <astro_ecfeci_sys.h>
#include <astro_keplerian.h>

#include <eomx.h>


/**
 * See eomx.h
 */
void eomx_print_orbits(const std::unordered_map<std::string,
                             std::shared_ptr<eom::Ephemeris>>& ephemerides,
                       const std::string& file_name,
                       const eom::EcfEciSys& f2iSys)
{
    // Print orbits to stdout in Cartesian and Keplerian
    // using computational reference frame
  std::cout << "\n\nGenerated Orbits";
  for (const auto& [name, eph] : ephemerides) {
    std::cout << "\n  " << name;
    std::cout << "\n  " << eph->getEpoch().to_string() << "    GCRF";
    eom::Keplerian oeCart(eph->getStateVector(eph->getEpoch(),
                                              eom::EphemFrame::eci));
    std::cout << oeCart;
  }

    // If a valid filename is supplied, write orbit definitions
    // using multiple useful conventions
  if (file_name.length() > 0) {
    std::ofstream fout(file_name);
    if (fout.is_open()) {
      for (const auto& [name, eph] : ephemerides) {
        fout << "\n\n  " << name << "  " << eph->getEpoch().to_string();
        fout << "\nITRF";
        Eigen::Matrix<double, 6, 1> itrf =
            eph->getStateVector(eph->getEpoch(), eom::EphemFrame::ecf);
          // Since ITRF is the best choice for transmission of state
          // vectors, use higher precision when writing to file
        fout << std::fixed <<
                std::setprecision(6) <<
                "\n    {" << phy_const::km_per_du*itrf(0) << "  " <<
                             phy_const::km_per_du*itrf(1) << "  " <<
                             phy_const::km_per_du*itrf(2) << "} km" <<
                std::setprecision(9) <<
                "\n    {" <<
                phy_const::km_per_du*itrf(3)*phy_const::tu_per_sec << "  " <<
                phy_const::km_per_du*itrf(4)*phy_const::tu_per_sec << "  " <<
                phy_const::km_per_du*itrf(5)*phy_const::tu_per_sec <<
                "} km/sec";

        fout << "\nGCRF";
        Eigen::Matrix<double, 6, 1> gcrf = 
            eph->getStateVector(eph->getEpoch(), eom::EphemFrame::eci);
        eom::Keplerian gcrfE(gcrf);
        fout << gcrfE;
        fout << "\nJ2000";
        Eigen::Matrix<double, 6, 1> j2000;
        j2000.block<3,1>(0,0) = f2iSys.gcrf2j2000(gcrf.block<3,1>(0,0));
        j2000.block<3,1>(3,0) = f2iSys.gcrf2j2000(gcrf.block<3,1>(3,0));
        eom::Keplerian j2000E(j2000);
        fout << j2000E;

        fout << "\nteme";
        Eigen::Matrix<double, 6, 1> teme =
            f2iSys.ecf2teme(eph->getEpoch(), itrf.block<3,1>(0,0),
                                             itrf.block<3,1>(3,0));
        eom::Keplerian temeE(teme);
        fout << temeE;
      }
    } else {
      std::cerr << "\n\n  Invalid Orbit Summary Filename:  " <<
                   file_name << '\n';
    }
  }
}

