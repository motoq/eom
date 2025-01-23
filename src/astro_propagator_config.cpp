/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_propagator_config.h>

#include <cal_julian_date.h>

namespace eom {

void PropagatorConfig::setStartStopTime(const JulianDate& jdStart,
                                        const JulianDate& jdStop)
{
  m_jdStart = jdStart;
  m_jdStop = jdStop;
}


void PropagatorConfig::setPropagator(Propagator integration_method)
{
  m_propagator = integration_method;
}


void PropagatorConfig::setGravityModel(GravityModel gravity_model)
{
  m_gravity_model = gravity_model;
}


void PropagatorConfig::setSunGravityModel(SunGravityModel sun_gravity)
{
  m_sun_gravity = sun_gravity;
}


void PropagatorConfig::setMoonGravityModel(MoonGravityModel moon_gravity)
{
  m_moon_gravity = moon_gravity;
}


void  PropagatorConfig::enableOtherGravityModels() noexcept
{
  m_other_gravity = true;
}


void  PropagatorConfig::disableOtherGravityModels() noexcept
{
  m_other_gravity = false;
}


void PropagatorConfig::setDegreeOrder(int degree, int order)
{
  m_degree = degree;
  m_order = order;
}


}
