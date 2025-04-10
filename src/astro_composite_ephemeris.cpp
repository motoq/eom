/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_composite_ephemeris.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <mth_index_mapper.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

CompositeEphemeris::CompositeEphemeris(
    const std::string& name,
    const std::vector<JulianDate>& handover_times,
    const std::vector<std::shared_ptr<eom::Ephemeris>>& ephems) :
                                       m_name {name},
                                       m_ephemerides {ephems}
{
  using namespace std::string_literals;
    // Ensure number of handover times corresponds to number
    // of ephemeris sources
  if (handover_times.size()  !=  ephems.size() - 1) {
    throw std::invalid_argument(
        "CompositeEphemeris::CompositeEphemeris():  " +
        std::to_string(handover_times.size()) +
        " handover times is not compatible with " +
        std::to_string(ephems.size()) + " ephemeris sources");
  }
  auto jdMin = m_ephemerides.front()->getBeginTime();
  auto jdMax = m_ephemerides.back()->getEndTime();
    // Ensure ephemeris sources cover handover times
  for (unsigned int ii=0; ii<handover_times.size(); ++ii) {
    if (handover_times[ii] <= jdMin  ||  jdMax <= handover_times[ii]) {
    throw std::invalid_argument(
        "CompositeEphemeris::CompositeEphemeris():  "s +
        "Handover time outside range of supplied ephemerides");
    }
  }
    // Create time blocks - IndexMapper with further validate
  std::vector<std::pair<JulianDate, JulianDate>> time_blocks;
  auto jd_share = jdMin;
  for (const auto& jd : handover_times) {
    time_blocks.emplace_back(jd_share, jd);
    jd_share = jd;
  }
  time_blocks.emplace_back(jd_share, jdMax);

  m_ndxr = std::make_unique<IndexMapper<JulianDate>>(time_blocks);
}


std::string CompositeEphemeris::getName() const
{
  return m_name;
}


JulianDate CompositeEphemeris::getEpoch() const
{
  return m_ephemerides.front()->getEpoch();
}


JulianDate CompositeEphemeris::getBeginTime() const
{
  return m_ephemerides.front()->getBeginTime();
}


JulianDate CompositeEphemeris::getEndTime() const
{
  return m_ephemerides.back()->getEndTime();
}


Eigen::Matrix<double, 6, 1>
CompositeEphemeris::getStateVector(const JulianDate& jd,
                                   EphemFrame frame) const 
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("CompositeEphemeris::getStateVector() - bad time");
  }
  return m_ephemerides[ndx]->getStateVector(jd, frame);
}


Eigen::Matrix<double, 3, 1>
CompositeEphemeris::getPosition(const JulianDate& jd,
                                EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("CompositeEphemeris::getPosition() - bad time");
  }
  return m_ephemerides[ndx]->getPosition(jd, frame);
}


}
