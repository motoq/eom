/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <stdexcept>
#include <string>

#include <cal_greg_date.h>

// Internal constants
namespace {
  static constexpr int gregyear {1582};        ///< Gregorian Calendar epoch
  static constexpr int maxyear  {9999};
  static constexpr int jan  {1};
  static constexpr int feb  {2};
  static constexpr int dec {12};
}

namespace eom {


GregDate::GregDate(int year, int month, int day)
{
  this->set(year, month, day);
}


GregDate::GregDate(std::string year, std::string month, std::string day)
{
  this->set(std::stoi(year), std::stoi(month), std::stoi(day));
}


void GregDate::set(std::string year, std::string month, std::string day)
{
  this->set(std::stoi(year), std::stoi(month), std::stoi(day));
}


void GregDate::set(int year, int month, int day)
{
    // Days in a standard (non leap year) month.
  static constexpr int nDaysInMonth[] = {  0,
                                          31, 28, 31, 30,
                                          31, 30, 31, 31,
                                          30, 31, 30, 31 };

    // Check year and month ranges
  if (year < gregyear  ||  year > maxyear) {
    std::string bad_year = "Invalid Year: " + std::to_string(year);
    throw std::invalid_argument(bad_year);
  } else if (month < jan  ||  month > dec) {
    std::string bad_month = "Invalid Month: " + std::to_string(month);
    throw std::invalid_argument(bad_month);
  } else {
    int offset = 0;
    if (month == feb  &&  isLeapYear(year)) {
      offset = 1;
    }

    int dim = nDaysInMonth[month] + offset;
    if (day < 1  ||  day > dim) {
      std::string bad_day = "Invalid Day: " + std::to_string(day);
      throw std::invalid_argument(bad_day);
    }
  }

  this->year = year;
  this->month = month;
  this->day = day;
}


bool GregDate::isLeapYear(int year) const noexcept
{
  bool leap = false;

  if (year < gregyear) {
    if (year%4 == 0) {
      leap = true;
    }
  } else {
    if (year%4 == 0) {
      if (year%100 == 0) {
        if (year%400 == 0) {
            leap = true;
        }
      } else {
        leap = true;
      }
    }
  }
  return leap;
}


}
