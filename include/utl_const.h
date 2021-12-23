/*
 * Copyright 2016 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_CONST_H
#define UTL_CONST_H

/**
 * Mathematical constants such as pi or conversion factors based on standards
 * (vs. emperical "constants" such as GM).
 *
 * @author  Kurt Motekew
 * @date    20160314
 */
namespace utl_const {

constexpr double pi   {3.141592653589793238462643383279502884};
constexpr double pio2 {pi/2.0};
constexpr double tpi  {2.0*pi};

constexpr double arcsec_per_rad {180.0*3600.0/pi};         //! Per Radian
constexpr double mas_per_rad {1000.0*arcsec_per_rad};      //! Milliarcsecond

}

#endif

