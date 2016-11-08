#include "differential.h"
#include <iostream>
#include <cmath>
using namespace std;
void Euler::setTimeStep(double h)
{
  _dt = h;
}
void Euler::calculateOneStep(SolarSystem *system)
{
  system->computeInteraction();
  for (CelestialBody &body : system->_celestialBodies)
  {
    body.position() += _dt*body.velocity();
    body.velocity() += _dt*body.interaction();
  }
  return;
}
void Verlet::setTimeStep(double h)
{
  _dt = h;
}
void Verlet::calculateOneStep(SolarSystem *system)
{
  double dt2 = 0.5*_dt;
  system -> computeInteraction();
  for(CelestialBody &body : system->_celestialBodies)
  {
    body.velocity() += dt2*body.interaction();
    body.position() += _dt*body.velocity();
  }

  system -> computeInteraction();
  for(CelestialBody &body : system->_celestialBodies)
  {
    body.velocity() += dt2*body.interaction();
  }
  return;
}
void VerletRelativistic::setTimeStep(double h)
{
  _dt = h;
  return;
}
void VerletRelativistic::calculateOneStep(SolarSystem *system)
{
  double dt2 = 0.5*_dt;
  system -> computeInteractionRelativistic();
  for(CelestialBody &body : system->_celestialBodies)
  {
    body.velocity() += dt2*body.interaction();
    body.position() += _dt*body.velocity();
  }

  system -> computeInteractionRelativistic();
  for(CelestialBody &body : system->_celestialBodies)
  {
    body.velocity() += dt2*body.interaction();
  }
  return;
}
