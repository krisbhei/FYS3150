#include "celestialbody.h"
#include <cmath>
#include <iostream>
using namespace std;
CelestialBody::CelestialBody() //Assuming to be the Sun
{
  _position = vec3(0,0,0);
  _velocity = vec3(0,0,0);
  _mass = 1;
  _interaction.zeros();
  _id = "-";
}
CelestialBody::CelestialBody(vec3 pos, vec3 vel, double mass, char const* id)
{
  _position = pos;
  _velocity = vel;
  _mass = mass;
  _interaction.zeros();
  _id = id;
}

double CelestialBody::getMass()
{
  return _mass;
}

vec3 &CelestialBody::interaction()
{
  return _interaction;
}

vec3 &CelestialBody::position()
{
  return _position;
}
vec3 &CelestialBody::velocity()
{
  return _velocity;
}
