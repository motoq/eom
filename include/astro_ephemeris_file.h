/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_EPHEMERIS_FILE_H
#define ASTRO_EPHEMERIS_FILE_H

#include <string>

namespace eom {

/**
 * Supported ephemeris file formats
 */
enum class EphFileFormat {
  sp3c                            ///< NGS Standard Product 3 format
};

/**
 * Supported interpolation methods
 */
enum class EphInterpType {
  hermite
};


/**
 * Holds parameters defining an ephemeris file and how it should be
 * processed.
 *
 * @author  Kurt Motekew
 * @date    2023/01/08
 */
class EphemerisFile {
public:
  /**
   * Create ephemeris file definition given a filename.
   *
   * @param  name        Unique ephemeris identifier
   * @param  eph_file    Name of file with ephemeris
   * @param  eph_format  Ephemeris format
   * @param  eph_interp  Interpolation method to be used
   */
  EphemerisFile(const std::string& name,
                const std::string& eph_file,
                EphFileFormat eph_format,
                EphInterpType eph_interp) : m_name {name},
                                            m_eph_file {eph_file},
                                            m_eph_format {eph_format},
                                            m_eph_interp {eph_interp}
  {
  }

  /*
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getName() const noexcept { return m_name; }

  /*
   * @return  Filename containing ephemeris
   */
  std::string getEphFileName() const noexcept { return m_eph_file; }

  /**
   * @return  Ephemeris format/type
   */
  EphFileFormat getEphFileFormat() const noexcept { return m_eph_format; }

  /**
   * @return  Interpolation method
   */
  EphInterpType getEphInterpMethod() const noexcept { return m_eph_interp; }


private:
  std::string m_name;
  std::string m_eph_file;
  EphFileFormat m_eph_format;
  EphInterpType m_eph_interp;
};


}

#endif
