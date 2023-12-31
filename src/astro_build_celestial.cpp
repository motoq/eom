/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <astro_build.h>

namespace eom {

std::vector<state_vector_rec>
build_celestial(const std::string& name_prefix,
                const JulianDate& startTime,
                const JulianDate& stopTime)
{
    // Read binary .emb file
  std::string fname = name_prefix + ".emb";
  std::ifstream ifs(fname, std::ios::binary );
  if (!ifs.is_open()) {
    throw std::runtime_error("build_celestial() Can't open " + fname);
  }

  double dt_days;
  ifs.read(reinterpret_cast<char*>(&dt_days), sizeof(double));
  double km_per_au;
  ifs.read(reinterpret_cast<char*>(&km_per_au), sizeof(double));
  const int rec_size {8};
  double sv_rec[rec_size];                  // jdhi, jdlo, x, y, z, dx, dy, dz
  int nrec {0};
  std::vector<state_vector_rec> sv_recs;
  JulianDate jd1 = startTime + -2.0*dt_days;
  JulianDate jd2 = stopTime  +  2.0*dt_days;
  bool covered_set {false};
  while (ifs.read(reinterpret_cast<char*>(sv_rec), rec_size*sizeof(double))) {
    JulianDate jd(sv_rec[0], sv_rec[1]);
    if (jd1 <= jd  &&  jd  <=  jd2) {
      Eigen::Matrix<double, 3, 1> pos = {sv_rec[2], sv_rec[3], sv_rec[4]};
      Eigen::Matrix<double, 3, 1> vel = {sv_rec[5], sv_rec[6], sv_rec[7]};
        // Stored in units of AU and days.  Leave in J2000 frame
      pos *= phy_const::du_per_km*km_per_au;
      vel *= phy_const::du_per_km*km_per_au*phy_const::day_per_tu;
      sv_recs.emplace_back(jd, pos, vel);
    }
    if (jd2 <= jd) {
      covered_set = true;
    }
  }
  ifs.close();

  if (sv_recs.size() < 2) {
    covered_set = false;
  }
  if (!covered_set) {
    throw std::runtime_error("build_celestial() Ephemeris not covered" +
                              fname);
  }
  return sv_recs;
}


}
