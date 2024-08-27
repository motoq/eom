/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_JULIAN_DATE_H
#define CAL_JULIAN_DATE_H

#include <string>

#include <cal_const.h>
#include <cal_greg_date.h>
#include <cal_duration.h>

namespace eom {

/**
 * A Julian Date class designed to preserve precision by splitting the date
 * into high and low values.
 *
 * @author  Kurt Motekew
 * @date    20160314
 * @date    20210916  Removed fractional second to simplify logic
 * @date    20210929  Added use of Duration for + and +=
 */
class JulianDate {
public:
  /**
   * Initialize with January 1, 2000
   */
  JulianDate() {}

  /**
   * Initialize using Julian Date components.
   *
   * @param  jdDays   Days portion of Julian Date, or full JD if other
   *                  parameters are zero or not included.  Allowing
   *                  implicit conversion seems appropriate for this
   *                  class.
   * @param  jdFrac   Additional days and fraction of a day to add to
   *                  the Julian Date.
   */
  explicit JulianDate(double jdDays, double jdFrac = 0.0) :
                      jdHi {jdDays}, jdLo {jdFrac}
  {
  }

  /**
   * Initialize with a Gregorian date and time of day.  The time parameters
   * are not bounded typical clock 12/24 and 60 base systems - as many
   * hours, minutes, and seconds from the gd, +/-, are allowed and
   * incorporated into the internal JD representation.
   *
   * @param  gd   Gregorian Date
   * @param  hr   Hours from gd
   * @param  min  Minutes from gd and hr
   * @param  sec  Seconds from gd, hr, and min
   */
  JulianDate(const GregDate& gd, int hr = 0, int min = 0, double sec = 0.0);

  /**
   * See the constructor with the same signature
   */
  void set(const GregDate& gd, int hr = 0, int min = 0, double sec = 0.0);

  /**
   * Returns the Julian date as a single double precision value.
   * For a single value representing time, the MJD preserves more
   * precision.
   *
   * @return   Julian Date, scalar form
   */
  double getJd() const noexcept
  {
    return jdLo + jdHi;
  }

  /**
   * @return   Modified Julian Date, scalar
   */
  double getMjd() const noexcept
  {
    return jdLo + (jdHi - cal_const::mjd);
  }

  void setMjd2000(double mjd2000) noexcept;

  /**
   * @return  Days since Jan 1, 2000
   */
  double getMjd2000() const noexcept
  {
    return jdLo + (jdHi - cal_const::j2000);
  }

  /**
   * @return  Julian Centuries
   */
  double getJulianCenturies() const noexcept
  {
    return (jdHi - 2451545.0 + jdLo)/36525.0;
  }

  /**
   * This method allows for compatibility with external libraries
   * that are designed to make use of a two part Julian date.
   *
   * @return   Large portion of the Julian Date, on the order of
   *           2,400,000.  Ideally, this would be the Julian Date
   *           corresponding to noon, but this is not necessary,
   *           depending on the value of jdLo.  Units of days.
   */
  double getJdHigh() const noexcept
  {
    return jdHi;
  }

  /**
   * This method allows for compatibility with external libraries
   * that are designed to make use of a two part Julian date.
   *
   * @return   Small portion of the Julian Date, typically on the order
   *           of a day or fraction of a day.  Units of days.
   */
  double getJdLow() const noexcept
  {
    return jdLo;
  }

  /**
   * Update this Julian Date by the given number of days.
   *
   * @param   days   Days to add to this Julian Date (or subtracted, if
   *                 negative).
   *
   * @return  This + days
   */
  JulianDate& operator+=(double days) noexcept;

  /**
   * Return a Julian date adjusted by the given number of days.
   *
   * @param   days   Days to add (or subtract, if negative).
   *
   * @return   Copy of this JulianDate + days
   */
  JulianDate operator+(double days) const noexcept;

  /**
   * Update this Julian Date by the given duration.
   *
   * @param   dur   Duration to add to this Julian Date (or subtracted, if
   *                negative).
   */
  JulianDate& operator+=(const Duration& dur) noexcept;

  /**
   * Return a Julian date adjusted by the given duration.
   *
   * @param   dur  Duration to add (or subtrack, if negative).
   *
   * @return   Copy of this JulianDate + the input Duration
   */
  JulianDate operator+(const Duration& dur) const noexcept;

  /**
   * @param  Julian date to subtract
   *
   * @return  The time difference, in days, between this JD and the
   *          input JD.  This - input.
   */
  double operator-(const JulianDate& jd) const noexcept
  {
    return jdHi - jd.jdHi + (jdLo - jd.jdLo);
  }

  /**
   * @param  Julian date to compare to
   *
   * @return  true if this JD is less than the other, jd1 < jd2.
   */
  bool operator<(const JulianDate& jd) const noexcept
  {
    return jdHi - jd.jdHi + (jdLo - jd.jdLo) < 0.0;
  }

  /**
   * @param  Julian date to compare to
   *
   * @return  true if this JD is less than or equal to the other, jd1 <= jd2.
   */
  bool operator<=(const JulianDate& jd) const noexcept
  {
    return jdHi - jd.jdHi + (jdLo - jd.jdLo) <= 0.0;
  }

  /**
   * @return   Gregorian Date and time as a string.  Time is
   *           in base 24:60:60 format:  yyyy/mm/dd hh:mm:ss.00
   */
  std::string to_str() const;

    /**
     * @return   Gregorian Date and time as a string.  Time is
     *           in base 24:60:60 format:  dd month yyyy hh:mm:ss.000000
     */
  std::string to_dmy_str() const;

private:
  void normalize();

  /*
   * Converts this JulianDate to a Gregorian date and time of day.
   *
   *   @param   year     Four digit year                                (output)
   *   @param   month    Month, 1-12                                    (output)
   *   @param   day      Day, 1-{28,29,30,31}                           (output)
   *   @param   hour     Hour of day, 0 <= hour < 24                    (output)
   *   @param   minutes  Minutes, 0 <= minutes < 60                     (output)
   *   @param   seconds  Seconds, 0 <= seconds < 60                     (output)
   */
  void jd2gd(int& year, int& month, int& day,
             int& hour, int& minutes, double& seconds,
             int iter = 0) const;

  double jdHi {cal_const::j2000};           // Days
  double jdLo {0.0};                        // Fraction of a day
};


}

#endif
