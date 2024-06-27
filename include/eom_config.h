/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_CONFIG_H
#define EOM_CONFIG_H

#include <deque>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <cal_duration.h>
#include <cal_julian_date.h>

namespace eom_app {


/**
 * A class containing EOM scenario configuration parameters and methods
 * for parsing them from tokenized string representations.
 *
 * @author  Kurt Motekew
 * @date    20210922
 */
class EomConfig {
public:
  /**
   * Initialize with default values.  The default settings are valid,
   * but not very useful.
   */
  EomConfig() {}

  /**
   * @param  tokens  Tokenized parameters representing the start time of
   *                 the simulation
   */
  void setStartTime(std::deque<std::string>& tokens);

  /**
   * @return  The start time of the simulation
   */
  eom::JulianDate getStartTime() const noexcept { return jdStart; }

  /**
   * @return  The end time of the simulation
   */
  eom::JulianDate getStopTime() const noexcept { return jdStop; }

  /**
   * @param  tokens  Tokenized parameters indicating the duration of the
   *                 scenario
   */
  void setDuration(std::deque<std::string>& tokens);

  /**
   * @param  tokens  Tokenized parameters indicating the number of assumed leap
   *                 seconds
   */
  void setLeapSeconds(std::deque<std::string>& tokens);

  /**
   * @param  tokens  Tokenized parameters representing the time update rate
   *                 that will be used when storing reduction parameters
   */
  void setEcfEciRate(std::deque<std::string>& tokens);

  /**
   * @return  The duration between stored reduction parameters
   */
  eom::Duration getEcfEciRate() const noexcept { return dtEcfEci; }

  /**
   * @param  tokens  Label used to determine and set the conversion
   *                 factor from input/output units to internal angle
   *                 units of radians.
   */
  void setIoPerRad(std::deque<std::string>& tokens);

  /**
   * @return  Get input/output units per radians.
   */
  double getIoPerRad() const noexcept { return io_per_rad; }

  /**
   * @return  Label describing output and expected input distance units
   */
  std::string getIoAngleUnits() const noexcept { return units_angle; }

  /**
   * @param  tokens  Label used to determine and set the conversion
   *                 factor from input/output units to internal distance
   *                 units.
   */
  void setIoPerDu(std::deque<std::string>& tokens);

  /**
   * @return  Get input/output units per internal distance units.
   */
  double getIoPerDu() const noexcept { return io_per_du; }

  /**
   * @return  Label describing output and expected input distance units
   */
  std::string getIoDistansUnits() const noexcept { return units_distance; }

  /**
   * @param  tokens  Label used to determine and set the conversion
   *                 factor from input/output units to internal time
   *                 units.
   */
  void setIoPerTu(std::deque<std::string>& tokens);

  /**
   * @return  Get input/output units per internal time units.
   */
  double getIoPerTu() const noexcept { return io_per_tu; }

  /**
   * @return  Label describing output and expected input time units
   */
  std::string getIoTimeUnits() const noexcept { return units_time; }

  /**
   * @param  tokens  Tokenized parameters indicating the duration of the
   *                 scenario
   */
  void setOutputRate(std::deque<std::string>& tokens);

  /**
   * @return  The period of time between outputs
   */
  eom::Duration getOutputRate() const noexcept { return dtOut; }

  /**
   * @param  orbit_name  Name of orbit to add to the list of orbits that
   *                     are to be created
   */
  void addPendingOrbit(const std::string& orbit_name);

  /**
   * @param  orbit_name  Name of orbit to look up
   *
   * @return  true if orbit_name has been entered as an orbit to be
   *          created
   */
  bool pendingOrbit(const std::string& orbit_name) const {
    return m_orbit_names.count(orbit_name) > 0;
  }


  /**
   * @param  name  Name of celestial object for which ephemeris should
   *               be loaded
   */
  void addCelestial(const std::string& name);

  /*
   * @return  Names of all celestial objects to load
   */
  std::vector<std::string> getCelestials() const;

  /**
   * @return  If an error was encountered while building the scenario,
   *          the return value is false.  Call getError().
   */
  bool isValid() const noexcept { return valid; }

  /**
   * @return  A string indicating an error state.  The string should be
   *          empty if no errors are present.
   */
  std::string getError() { return error_string; }

private:
  std::string units_angle {"Radians"};
  std::string units_distance {"DU"};
  std::string units_time {"TU"};
  double io_per_rad {1.0};
  double io_per_du {1.0};
  double io_per_tu {1.0};
  bool valid {true};
  bool epoch_set {false};
  bool f2i_rate_set {false};
  bool leapsec_set {false};
  std::string error_string {""};
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  eom::Duration dtEcfEci;
  eom::Duration dtOut;

  std::set<std::string> m_orbit_names;
  std::vector<std::string> m_celestial_names;
  
};


/**
 * Override output stream as non-member function.
 */
std::ostream& operator<<(std::ostream& out, const EomConfig& cfg);

}

#endif
