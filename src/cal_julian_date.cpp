/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cal_julian_date.h>

#include <string>
#include <cmath>

#include <cal_const.h>
#include <cal_greg_date.h>
#include <cal_duration.h>

static double gd2jd(int iyear, int imonth, int iday);

namespace eom {

JulianDate::JulianDate(const GregDate& gd, int hr, int min, double sec)
{
  this->set(gd, hr, min, sec);
}


void JulianDate::set(const GregDate& gd, int hr, int min, double sec)
{
  jdHi = gd2jd(gd.getYear(), gd.getMonth(), gd.getDay());
  jdLo = cal_const::day_per_min*(60*hr + min + cal_const::min_per_sec*sec);

  if (jdLo != 0.0) {
    double more_days = static_cast<long>(jdLo);
    jdHi += more_days;
    jdLo -= more_days;
  }
}


void JulianDate::setMjd2000(double mjd2000) noexcept
{
  jdHi = cal_const::j2000;
  jdLo = mjd2000;
}


/*
 * Break into low and high for more precision
 */
JulianDate& JulianDate::operator+=(double days) noexcept
{
  //double full_days = static_cast<long>(days);

  //jdHi += full_days;
  //jdLo += (days - full_days);
  jdLo += days;

  return *this;
}


JulianDate JulianDate::operator+(double days) const noexcept
{
  JulianDate jd {*this};
  jd += days;

  return jd;
}


/*
 * Break into low and high for more precision
 */
JulianDate& JulianDate::operator+=(const Duration& dur) noexcept
{
  double days {dur.getDays()};
  //double full_days = static_cast<long>(days);

  //jdHi += full_days;
  //jdLo += (days - full_days);
  jdLo += days;

  return *this;
}


JulianDate JulianDate::operator+(const Duration& dur) const noexcept
{
  JulianDate jd {*this};
  double days {dur.getDays()};
  jd += days;

  return jd;
}


/**
 * Note: Seconds are truncated to 1/100 th of a second
 */
std::string JulianDate::to_str() const
{
  int year, month, day, hour, minutes;
  double seconds; 
  jd2gd(year, month, day, hour, minutes, seconds);

    // Stop *printf from rounding to silly values like xx:xx:60 sec
  seconds = 0.01*static_cast<long>(100.0*seconds);

  char buf[32];
  snprintf(buf, sizeof(buf), "%4i/%02i/%02i %02i:%02i:%05.2f",
                              year, month, day, hour, minutes, seconds);

  std::string dts{buf};
  return dts;
}


/**
 * Note: Seconds are truncated to 1/100 th of a second
 */
std::string JulianDate::to_dmy_str() const
{
  int year, month, day, hour, minutes;
  double seconds; 
  jd2gd(year, month, day, hour, minutes, seconds);
  GregDate gd(year, month, day);

    // Stop *printf from rounding to silly values like xx:xx:60 sec
  seconds = 0.000001*static_cast<long>(1000000.0*seconds);

  char buf[32];
  snprintf(buf, sizeof(buf), "%i %s %i %02i:%02i:%09.6f", 
                              day, gd.getMonthStr().c_str(), year,
                              hour, minutes, seconds);

  std::string dts{buf};
  return dts;
}


/*
 * Sets internal jdHi to be at the beginning of the day (ends in .5) with
 * seconds rolled into jdLo and jdLo positive and less than one day.
 */
void JulianDate::normalize()
{
    // First move jdHi to noon and move seconds to jdLo
    // (floor goes back to the previous day, then add 1 for today at noon)
  double tmp {std::floor(jdHi) + 1.0};
  jdLo += (jdHi - tmp);
  jdHi = tmp;

    // Now, push whole days of jdLo into jdHi - jdLo sign may be +/-
  double more_days = static_cast<long>(jdLo);
  jdHi  += more_days;
  jdLo -= more_days;

    // Move jdHi to 00:00:00
  jdHi  -= 0.5;
  jdLo += 0.5;

    // Keep jdLo positive and in the range 00:00:00 <= jdLo < 24:00:00
  if (jdLo < 0.0) {
    jdLo += 1.0;
    jdHi  -= 1.0;
  } else if (jdLo > 1.0) {
    jdLo -= 1.0;
    jdHi  += 1.0;
  }
}


/*
 * Converts this JulianDate to a Gregorian date and time of day.
 * This algorithm comes from the U.S. Naval Observatory website, the
 * astronomical Applications Department.  They got the formula from the
 * Fliegel reference.
 *
 * Note: The USNO website had this formula back around the year 2000.
 * As the site is currently undergoing "modernization", most likely by
 * useless agile deveolpers, everything that made the site great before
 * will most likely be lost forever.  The interested reader is encouraged
 * to pick up a copy of the Explanatory Supplement to the Astronomical
 * Almanac, if interested in the source of these formulas.  Also, the
 * IAU publishes easy to use C and FORTRAN code that performs these
 * calculations.
 *
 *   @param   year     Four digit year                                 (output)
 *   @param   month    Month, 1-12                                     (output)
 *   @param   day      Day, 1-{28,29,30,31}                            (output)
 *   @param   hour     Hour of day, 0 <= hour < 24                     (output)
 *   @param   minutes  Minutes, 0 <= minutes < 60                      (output)
 *   @param   seconds  Seconds, 0 <= seconds <= 60 (60th for leapsec)  (output)
 */
void JulianDate::jd2gd(int& year, int& month, int& day,
                       int& hour, int& minutes, double& seconds) const
{
    // Must get jdHi and jdLo in proper form first
  JulianDate tmpJd = *this;
  tmpJd.normalize();

  long jd {1L + static_cast<long>(tmpJd.jdHi)};
  long i, j, k, m, n;

  m = jd + 68569L;
  n = 4L*m/146097L;
  m = m - (146097L*n + 3L)/4L;
  i = 4000L*(m + 1L)/1461001L;
  m = m - 1461L*i/4L + 31L;
  j = 80L*m/2447L;
  k = m - 2447L*j/80L;            // day
  m = j/11L;
  j = j + 2L - 12L*m;             // month
  i = 100L*(n - 49L) + i + m;     // year
    //
  year  = i;
  month = j;
  day   = k;

  double hours_left {cal_const::hr_per_day*tmpJd.jdLo};
  hour = static_cast<int>(hours_left);

  double minutes_left {60.0*(hours_left - hour)};
  minutes = static_cast<int>(minutes_left);

  seconds = 60.0*(minutes_left - minutes);
}


}


/*
 * Calculates the Julan Date jd from the year, month, and day values.
 * This formula comes from the U.S. Naval Observatory website, the
 * astronomical Applications Department.  They got the formula from the
 * Fliegel reference.  It should be valid for any  dates that result in
 * a Julian date greater than zero.
 *
 * Note: The USNO website had this formula back around the year 2000.
 * As the site is currently undergoing "modernization", most likely by
 * useless agile deveolpers, everything that made the site great before
 * will most likely be lost forever.  The interested reader is encouraged
 * to pick up a copy of the Explanatory Supplement to the Astronomical
 * Almanac, if interested in the source of these formulas.  Also, the
 * IAU publishes easy to use C and FORTRAN code that performs these
 * calculations.
 *
 * @return   The JD at 00:00:00 of the input Gregorian Date
 */
static double gd2jd(int iyear, int imonth, int iday)
{
  long year {iyear};
  long month {imonth};
  long day {iday};

  long jd {day - 32075L +
                 1461L*(year + 4800L + (month - 14L)/12L)/4L +
                  367L*(month - 2L - (month-14L)/12L*12L)/12L -
                    3L*((year + 4900L + (month -14L)/12L)/100L)/4L};

  return static_cast<double>(jd) - 0.5;
}
