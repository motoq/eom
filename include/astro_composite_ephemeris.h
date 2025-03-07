/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_COMPOSITE_EPHEMERIS_H
#define ASTRO_COMPOSITE_EPHEMERIS_H

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <mth_index_mapper.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>


namespace eom {

/**
 * Ephemeris composed of multiple Ephemeris objects stitched together.
 * Selection - switchover from one Ephemeris object to another - is
 * based on time.  Continuity is not guaranteed/checked by this class.
 * If continuity is desired, it must be ensured while populating this
 * object.
 *
 * @author Kurt Motekew  2025/02/15
 */
class CompositeEphemeris : public Ephemeris {
public:
  /**
   * Initializes given multiple ephemeris sources and a schedule
   * for which to transition from one to the other.
   *
   * @param  name            Unique ephemeris identifier
   * @param  handover_times  Handover schedule from one ephemeris source
   *                         to the next.  There must be n-1 time for n
   *                         ephemerides supplied.  The times must be
   *                         sequential and covered by the ephemerides.
   *                         For example, if the first handoff occurs
   *                         at time t, then the first ephemeris source
   *                         must have an end time covering t while
   *                         the second ephemeris source must have a
   *                         begin time that occurs at or before t.
   * @param  ephems          Ephemeris sources for which the shared pointers
   *                         will be copied into this class for use.
   */
  CompositeEphemeris(
      const std::string& name,
      const std::vector<JulianDate>& handover_times,
      const std::vector<std::shared_ptr<eom::Ephemeris>>& ephems);

  /**
   * @return  Ephemeris identifier
   */
  std::string getName() const override;

  /**
   * @return  Epoch for first set of ephemerides defining this composite
   *          set.
   */
  JulianDate getEpoch() const override;

  /**
   * @return  Earliest time for which ephemeris can be extracted
   */
  JulianDate getBeginTime() const override;

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override;

  /**
   * @param  jd     Time for which to return a state vector, UTC
   * @param  frame  Reference frame of returned state vector
   *
   * @return  Cartesian position and velocity state vector, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame) const override;

  /**
   * @param  jd     Time for which to return a position vector, UTC
   * @param  frame  Reference frame of returned position vector
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate&,
                                          EphemFrame) const override;

private:
  std::string m_name;
  std::vector<std::shared_ptr<eom::Ephemeris>> m_ephemerides;
  std::unique_ptr<IndexMapper<JulianDate>> m_ndxr {nullptr};
};


}

#endif
