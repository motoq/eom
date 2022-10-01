/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cal_greg_date.h>

#include <string>
#include <stdexcept>

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
    throw std::invalid_argument("GregDate::set() Invalid Year: " +
                                std::to_string(year));
  } else if (month < jan  ||  month > dec) {
    throw std::invalid_argument("GregDate::set() Invalid Month: " +
                                std::to_string(month));
  } else {
    int offset = 0;
    if (month == feb  &&  isLeapYear(year)) {
      offset = 1;
    }

    int dim = nDaysInMonth[month] + offset;
    if (day < 1  ||  day > dim) {
      throw std::invalid_argument("GregDate::set() Invalid Day: " +
                                  std::to_string(day));
    }
  }

  this->year = year;
  this->month = month;
  this->day = day;
}


std::string GregDate::getMonthStr() const noexcept
{
  switch (month) {
    case 1:  return "Jan";
    case 2:  return "Feb";
    case 3:  return "Mar";
    case 4:  return "Apr";
    case 5:  return "May";
    case 6:  return "Jun";
    case 7:  return "Jul";
    case 8:  return "Aug";
    case 9:  return "Sep";
    case 10: return "Oct";
    case 11: return "Nov";
    case 12: return "Dec";
    default: return "";
  }
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
