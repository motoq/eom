/*
 * Copyright 2022, 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_propagator_config.h>

#include <cal_duration.h>
#include <cal_julian_date.h>

namespace eom {

PropagatorConfig::PropagatorConfig()
{
}


PropagatorConfig::PropagatorConfig(PropagatorType prop_type) :
                                   m_prop_type {prop_type}
{
}


PropagatorType PropagatorConfig::getPropagatorType() const noexcept
{
  return m_prop_type;
}


void PropagatorConfig::setStartStopTime(const JulianDate& jdStart,
                                        const JulianDate& jdStop)
{
  m_jdStart = jdStart;
  m_jdStop = jdStop;
}


JulianDate PropagatorConfig::getStartTime() const noexcept
{
  return m_jdStart;
}


JulianDate PropagatorConfig::getStopTime() const noexcept
{
  return m_jdStop;
}


void PropagatorConfig::setPropagator(Propagator integration_method)
{
  m_propagator = integration_method;
}


Propagator PropagatorConfig::getPropagator() const noexcept
{
  return  m_propagator;
}


void PropagatorConfig::setStepSize(const Duration& dt)
{
  m_dt = dt;
}


Duration PropagatorConfig::getStepSize() const noexcept
{
  return  m_dt;
}


void PropagatorConfig::setGravityModel(GravityModel gravity_model)
{
  m_gravity_model = gravity_model;
}


GravityModel PropagatorConfig::getGravityModel() const noexcept
{
  return m_gravity_model;
}


void PropagatorConfig::setSunGravityModel(SunGravityModel sun_gravity)
{
  m_sun_gravity = sun_gravity;
}


SunGravityModel PropagatorConfig::getSunGravityModel() const noexcept
{
  return m_sun_gravity;
}


void PropagatorConfig::setMoonGravityModel(MoonGravityModel moon_gravity)
{
  m_moon_gravity = moon_gravity;
}


MoonGravityModel PropagatorConfig::getMoonGravityModel() const noexcept
{
  return m_moon_gravity;
}


void  PropagatorConfig::enableOtherGravityModels() noexcept
{
  m_other_gravity = true;
}


void  PropagatorConfig::disableOtherGravityModels() noexcept
{
  m_other_gravity = false;
}


bool PropagatorConfig::otherGravityModelsEnabled() const noexcept
{
  return m_other_gravity;
}


void PropagatorConfig::setDegreeOrder(int degree, int order)
{
  m_degree = degree;
  m_order = order;
}


int PropagatorConfig::getDegree() const noexcept
{
  return m_degree;
}


int PropagatorConfig::getOrder() const noexcept
{
  return m_order;
}



void PropagatorConfig::setSrpModel(SrpModel srp_model)
{
  m_srp_model = srp_model;
}


SrpModel PropagatorConfig::getSrpModel()
{
  return m_srp_model;
}


void PropagatorConfig::setReflectivity(double cr)
{
  m_cr = cr;
}


double PropagatorConfig::getReflectivity()
{
  return m_cr;
}


void PropagatorConfig::setAreaOverMass(double aom)
{
  m_aom = aom;
}


double PropagatorConfig::getAreaOverMass()
{
  return m_aom;
}


}
