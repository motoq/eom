/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_GREG_DATE_H
#define CAL_GREG_DATE_H

#include <string>

namespace eom {


/**
 * A simple Gregorian Date class primarily used to error verify that a
 * year, month, and day are valid.
 *
 * @author  Kurt Motekew
 * @date    20160314
 */
class GregDate {
public:
  /**
   * Initialized to the date of the Sputnik launch.
   */
  GregDate() {}

  /**
   * Initialize with integer representations of calendar values.
   *
   * @param  year   Four digit representation of the year
   * @param  month  1 <= month <= 12
   * @param  day    1 <= day <= {28, 29, 30, 31}, depending on the
   *                month/year combination.
   *
   * @throws invalid_argument
   */
  GregDate(int year, int month, int day);
 
  /**
   * Initialize with a string representations of calendar values.
   *
   * @param  year   Four digit representation of the year
   * @param  month  1 <= month <= 12
   * @param  day    1 <= day <= {28, 29, 30, 31}, depending on the
   *                 month/year combination.
   *
   * @throws invalid_argument
   */
  GregDate(std::string year, std::string month, std::string day);

  /**
   * Set with a string representations of calendar values.
   *
   * @param  year   Four digit representation of the year
   * @param  month  1 <= month <= 12
   * @param  day    1 <= day <= {28, 29, 30, 31}, depending on the
   *                 month/year combination.
   *
   * @throws invalid_argument
   */
  void set(std::string year, std::string month, std::string day);

  /**
   * Set date with integer representations of calendar values.
   *
   * @param  year   Four digit representation of the year
   * @param  month  1 <= month <= 12
   * @param  day    1 <= day <= {28, 29, 30, 31}, depending on the
   *                 month/year combination.
   *
   * @throws invalid_argument
   */
  void set(int year, int month, int day);

  /**
   * @return  Four digit representation of the year
   */
  double getYear() const noexcept { return year; }

  /**
   * @return  Month number, 1->12
   */
  double getMonth() const noexcept { return month; }

  /**
   * @return  Day of the month, 1->{28,29,30,31}
   */
  double getDay() const noexcept { return day; }

  /**
   * @return  Three letter text abbreviation of the month
   */
  std::string getMonthStr() const;

  /**
   * Determins if the given year is a leapyear.  The method used to determine
   * this is as follows.  If the year is divisible by 4 but not divisible by
   * 100 then the year is a leap year.  If the year is divisible by 4, and
   * divisible by 100, then it is a leap year only if it is also divisible by 
   * 400.  This method should be valid from the year 1582 forward since that is
   * when Pope Gregory XIII decided that this would be the method used to 
   * determie leap years as opposed to the previous method that only required 
   * the year be divisible by 4.  If the year is before 1582 then the year will
   * be considered a leap year if it is divisible by 4.  Not sure how far back 
   * that will work, but hey, nothing important happended back then.
   *
   * @param  year  Four digit representation of the year
   *
   * @return  If true, then this is a leap year
   */
  bool isLeapYear(int year) const noexcept;

private:
  int year {1957};
  int month {10};
  int day {4};
};


}

#endif
