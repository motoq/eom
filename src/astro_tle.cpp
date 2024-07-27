/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

#include <cal_greg_date.h>
#include <astro_tle.h>

#include <iostream>

/*
 * Given indexes (not offsets) into a string, pull the substring after
 * converting to an offset and length, and covert to a double.
 *
 * @param  ndx1     Index into TLE line of first character
 * @param  ndx2     Index into TLE line of last character
 * @param  tle_str  String of max length ndx2-1
 */
static double get_tle_double(int ndx1, int ndx2, const std::string& tle_str);

namespace eom {

Tle::Tle(const std::string& tle1, const std::string& tle2)
{
  if (tle1.length() < 63) {
    throw std::invalid_argument("Tle::Tle(): Invalid TLE line 1 length: " +
                                 std::to_string(tle1.length()));
  } else if (tle2.length() < 63) {
    throw std::invalid_argument("Tle::Tle(): Invalid TLE line 2 length: " +
                                 std::to_string(tle2.length()));
  }

  m_tle1 = tle1;
  m_tle2 = tle2;

  /*
   * ndx1 and ndx2 are TLE indexs for values.  This is being done to
   * make it easier to see how the data is being parsed given
   * documentation on the TLE format.
   */

  int ndx1, ndx2, offset, len;
    // Satellite (alphanumeric) number
  ndx1 = 3;
  ndx2 = 7;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  m_satn = tle1.substr(offset, len);
    // Satellite year
  ndx1 = 19;
  ndx2 = 20;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  auto oe_str = tle1.substr(offset, len);
  try {
    m_epoch_yr = std::stoi(oe_str);
      // Fucking Y2K nonsense 1/4 century late!? - seriously USSF!
      // Everyone has to change "Satellite Number" to be "Satellite Name".
      // Seems like a primo time to update the TLE format.
    GregDate gd;
    if (m_epoch_yr < 0  ||  m_epoch_yr > 99) {
      throw std::invalid_argument("");
    } else if (m_epoch_yr < 57) {
      gd.set(2000 + m_epoch_yr, 1, 1);
    } else {
      gd.set(1900 + m_epoch_yr, 1, 1);
    }
    m_epoch.set(gd);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Invalid Epoch Year: " +
                                 oe_str);
  }
    // Satellite day of year
  try {
    m_epoch_doy = get_tle_double(21, 32, tle1);
    if (m_epoch_doy < 1.0  ||  m_epoch_doy > 366.0) {
      throw std::invalid_argument("");
    }
    m_epoch += (m_epoch_doy - 1);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Invalid Epoch Day of Year: " +
                                 tle1);
  }
    // First time derivative of the mean motion
  try { 
    m_ndot = get_tle_double(34, 43, tle1);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): d(no)/dt: " + tle1);
  }
    // Second time derivative of the mean motion
  ndx1 = 45;
  ndx2 = 52;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  oe_str = tle1.substr(offset, len);
  std::string sign_str = oe_str.substr(0, 1);
  std::string decimal_str = ".";
  std::string mtsa_str = oe_str.substr(1, 5);
  std::string e_str = "e";
  std::string exp_str = oe_str.substr(6, 2);
  oe_str = sign_str + decimal_str + mtsa_str + e_str + exp_str;
  try { 
    m_nddot = std::stod(oe_str);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): d2(no)/dtdt: " + oe_str);
  }
    // B*
  ndx1 = 54;
  ndx2 = 61;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  oe_str = tle1.substr(offset, len);
  sign_str = oe_str.substr(0, 1);
  mtsa_str = oe_str.substr(1, 5);
  exp_str = oe_str.substr(6, 2);
  oe_str = sign_str + decimal_str + mtsa_str + e_str + exp_str;
  try { 
    m_bstar = std::stod(oe_str);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): B*: " + oe_str);
  }
    // Eph type
  ndx1 = 63;
  ndx2 = 63;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  oe_str = tle1.substr(offset, len);
  try { 
    m_eph_type = std::stoi(oe_str);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Invalid Eph Type: " +
                                 oe_str);
  }
    // Inclination
  try { 
    m_inclo = get_tle_double(9, 16, tle2);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Inclination: " + tle2);
  }
    // RAAN
  try { 
    m_nodeo = get_tle_double(18, 25, tle2);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): RAAN: " + tle2);
  }
    // Eccentricity
  ndx1 = 27;
  ndx2 = 33;
  offset = ndx1 - 1;
  len = ndx2 - offset;
  oe_str = tle2.substr(offset, len);
  oe_str = decimal_str + oe_str;
  try { 
    m_ecco = std::stod(oe_str);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): B*: " + oe_str);
  }
    // Mean anomaly
  try { 
    m_m0 = get_tle_double(44, 51, tle2);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Mean Anomaly: " + tle2);
  }
    // Mean motion
  try { 
    m_n0 = get_tle_double(53, 63, tle2);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Tle::Tle(): Mean Motion: " + tle2);
  }

}


std::string Tle::getLineOne() const
{
  return m_tle1;
}


std::string Tle::getLineTwo() const
{
  return m_tle2;
}


std::string Tle::getSatName() const
{
  return m_satn;
}


JulianDate Tle::getEpoch() const
{
  return m_epoch;
}


int Tle::getEpochYear() const
{
  return m_epoch_yr;
}


double Tle::getEpochDayOfYear() const
{
  return m_epoch_doy;
}


double Tle::getMeanMotionRate() const
{
  return m_ndot;
}


double Tle::getMeanMotionSecondRate() const
{
  return m_nddot;
}


double Tle::getBstar() const
{
  return m_bstar;
}


int Tle::getElementType() const
{
  return m_eph_type;
}


double Tle::getInclination() const
{
  return m_inclo;
}


double Tle::getRaan() const
{
  return m_nodeo;
}


double Tle::getEccentricity() const
{
  return m_ecco;
}


double Tle::getMeanAnomaly() const
{
  return m_m0;
}


double Tle::getMeanMotion() const
{
  return m_n0;
}


std::ostream& operator<<(std::ostream& out, const Tle& tle)
{
    // DOY already sanity checked - simple Y2K conversion
  auto year = 1900 + tle.getEpochYear();
  if (year < 1957) {
    year += 100;
  }
  return out << tle.getEpoch().to_str() << '\n' <<
                tle.getLineOne() << '\n' << tle.getLineTwo() << '\n' <<
                "  Designator:   " << tle.getSatName() << '\n' <<
                "  Year:DOY:     " << year << ':' <<
                                      std::fixed <<
                                      std::setprecision(8) <<
                                      tle.getEpochDayOfYear() << '\n' <<
                std::scientific <<
                std::setprecision(7) <<
                "  ndot/2:       " << tle.getMeanMotionRate() <<
                                      " rev/day^2\n" <<
                std::setprecision(4) <<
                "  nddot/6:      " << tle.getMeanMotionSecondRate() <<
                                 "  rev/day^3\n" <<
                "  B*:           " << tle.getBstar() << " 1/ER\n" <<
                "  Element Type: " << tle.getElementType() << '\n' <<
                std::fixed <<
                "  Inclination:  " << tle.getInclination() << "\u00B0\n" <<
                "  RAAN:         " << tle.getRaan() << "\u00B0\n" <<
                "  Eccentricity: " << tle.getEccentricity() << "\u00B0\n" <<
                "  Mean Anomaly: " << tle.getMeanAnomaly() << "\u00B0\n" <<
                "  Mean Motion:  " << tle.getMeanMotion() << "\u00B0\n";
}


}


static double get_tle_double(int ndx1, int ndx2, const std::string& tle_str)
{
  auto offset = ndx1 - 1;
  auto len = ndx2 - offset;
  auto oe_str = tle_str.substr(offset, len);

  return std::stod(oe_str);
}
