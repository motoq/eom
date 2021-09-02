/*
 * Copyright 2016 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <cmath>

#include <utl_greg_date.h>
#include <astro_julian_date.h>

/*
 * Set whole and fractional values - seconds will be lumped with
 * fractional value (converted to days).
 */
JulianDate::JulianDate(const JulianDate& jd)
{
  this->set(jd);
}


JulianDate::JulianDate(double jdDays, double jdFrac, double jdSeconds) :
                         jdHi{jdDays}, jdLow{jdFrac}, jdSec{jdSeconds}
{
}


JulianDate::JulianDate(const GregDate& gd, int hr, int min, double sec)
{
  this->set(gd, hr, min, sec);
}


/*
 * Set whole and fractional values - seconds will be lumped with
 * fractional value (converted to days).
 */
void JulianDate::set(const JulianDate& jd)
{
  jdHi = jd.jdHiVal();
  jdLow = jd.jdLowVal();
  jdSec = 0.0;
}


void JulianDate::set(double jdDays, double jdFrac, double jdSeconds)
{
  jdHi = jdDays;
  jdLow = jdFrac;
  jdSec = jdSeconds;
}


void JulianDate::set(const GregDate& gd, int hr, int min, double sec)
{
  jdSec = sec;
  jdHi = gd2jd(gd.year(), gd.month(), gd.day());
  jdLow = DAY_PER_MIN*(60*hr + min);

  if (jdLow != 0.0) {
    double more_days = static_cast<int>(jdLow);
    jdHi += more_days;
    jdLow -= more_days;
  }
}


double JulianDate::getJd() const
{
  if (jdSec == 0.0) {
    return jdLow + jdHi;
  } else {
    return DAY_PER_SEC*jdSec + jdLow + jdHi;
  }
}


double JulianDate::getMjd() const
{
  if (jdSec == 0.0) {
    return DAY_PER_MIN*jdLow + (jdHi - MJD);
  } else {
    return DAY_PER_SEC*jdSec + DAY_PER_MIN*jdLow + (jdHi - MJD);
  }
}


double JulianDate::getJdLow() const
{
  if (jdSec == 0.0) {
    return jdLow;
  } else {
    return DAY_PER_SEC*jdSec + jdLow;
  }
}


/*
 * Break into low and high for more precision
 */
JulianDate& JulianDate::operator+=(double days)
{
  int full_days = static_cast<int>(days);

  jdHi += full_days;
  jdLow += (days - full_days);

  return *this;
}


JulianDate JulianDate::operator+(double days)
{
  JulianDate jd {*this};
  jd += days;

  return jd;
}


double JulianDate::operator-(const JulianDate& jd)
{
  return this->jdHiVal() - jd.jdHiVal() + (this->jdLowVal() - jd.jdLowVal());
}


/*
std::ostream& operator<<(std::ostream& out, const JulianDate& jd)
{
  JulianDate jd2 {jd};
  out << jd2.to_str();
  return out;
}
*/


/**
 * Note: Seconds are truncated to 1/100 th of a second
 */
std::string JulianDate::to_str()
{
  int year, month, day, hour, minutes;
  double seconds; 
  jd2gd(year, month, day, hour, minutes, seconds);

    // Stop *printf from rounding to silly values like xx:xx:60 sec
  int isec = static_cast<int>(100.0*seconds);
  seconds = 0.01*static_cast<double>(isec);

  char buf[32];
  snprintf(buf, sizeof(buf), "%4i/%02i/%02i %02i:%02i:%05.2f",
                              year, month, day, hour, minutes, seconds);

  std::string dts{buf};
  return dts;
}


/*
 * @return   The JD at 00:00:00 of the input Gregorian Date
 */
double JulianDate::gd2jd(int year, int month, int day)
{
  int jd = day - 32075
               + 1461*(year+4800+(month-14)/12)/4+367*(month-2-(month-14)/12*12)
               / 12-3*((year+4900+(month-14)/12)/100)/4;

  return jd - 0.5;
}


/*
 * Sets internal jdHi to be at the beginning of the day (ends in .5) with
 * seconds rolled into jdLow and jdLow positive and less than one day.
 */
void JulianDate::normalize()
{
    // First move jdHi to noon and move seconds to jdLow
    // (floor goes back to the previous day, then add 1 for today at noon)
  double tmp = std::floor(jdHi) + 1.0;
  jdLow += (jdHi - tmp);
  if (jdSec != 0.0) {
    jdLow += DAY_PER_SEC*jdSec;
    jdSec = 0.0;
  }
  jdHi = tmp;

    // Now, push whole days of jdLow into jdHi - jdLow sign may be +/-
  double more_days = static_cast<int>(jdLow);
  jdHi  += more_days;
  jdLow -= more_days;

    // Move jdHi to 00:00:00
  jdHi  -= 0.5;
  jdLow += 0.5;

    // Keep jdLow positive and in the range 00:00:00 <= jdLow < 24:00:00
  if (jdLow < 0.0) {
    jdLow += 1.0;
    jdHi  -= 1.0;
  } else if (jdLow > 1.0) {
    jdLow -= 1.0;
    jdHi  += 1.0;
  }
}


/*
 * Converts this JulianDate to a Gregorian date and time of day.
 *
 *   @param   year     Four digit year                                 (output)
 *   @param   month    Month, 1-12                                     (output)
 *   @param   day      Day, 1-{28,29,30,31}                            (output)
 *   @param   hour     Hour of day, 0 <= hour < 24                     (output)
 *   @param   minutes  Minutes, 0 <= minutes < 60                      (output)
 *   @param   seconds  Seconds, 0 <= seconds <= 60 (60th for leapsec)  (output)
 */
void JulianDate::jd2gd(int& year, int& month, int& day,
                       int& hour, int& minutes, double& seconds)
{
    // Must get jdHi and jdLow in proper form first
  normalize();

  int jd = 1 + static_cast<int>(jdHi);
  int i, j, k, m, n;

  m = jd+68569;
  n = 4*m/146097;
  m = m-(146097*n+3)/4;
  i = 4000*(m+1)/1461001;
  m = m-1461*i/4+31;
  j = 80*m/2447;
  k = m-2447*j/80;             // day
  m = j/11;
  j = j+2-12*m;                // month
  i = 100*(n-49)+i+m;          // year
    //
  year  = i;
  month = j;
  day   = k;

  double hours_left = HR_PER_DAY*jdLow;
  hour = static_cast<int>(hours_left);

  double minutes_left = 60.0*(hours_left - hour);
  minutes = static_cast<int>(minutes_left);

  seconds = 60.0*(minutes_left - minutes);
}
