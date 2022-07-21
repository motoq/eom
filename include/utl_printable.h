/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_PRINTABLE_H
#define UTL_PRINTABLE_H

#include <ostream>

namespace eom {

/**
 * Interface for an EOM object that can print its contents
 *
 * @author  Kurt Motekew
 * @date    20220422
 */
class Printable {
public:
  virtual ~Printable() = default;
  Printable() = default;
  Printable(const Printable&) = default;              // vs. delete
  Printable& operator=(const Printable&) = delete;
  Printable(Printable&&) = delete;
  Printable& operator=(Printable&&) = delete;

  /**
   * Print contents to input stream
   */
  virtual void print(std::ostream& stream) const = 0;

};


}

#endif
