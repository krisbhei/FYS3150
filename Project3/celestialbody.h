#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H
#include "vec3.h"

class CelestialBody
{
public:

  /*
  CelestialBody()

  Constructur assumed to use for the Sun
  Sets _position,_velocity and _interaction to be zero.
  _mass is equal to one, and _id to an identifier such that
  the Sun is not plotted with zero velocity
  */
  CelestialBody();

  /*
  CelestialBody(vec3 pos, vec3 vel, double mass, char const* id)

  Input:
    - vec3 pos       : position of the body, assigned to _position
    - vec3 vel       : velocity of the body, assigned to _velocity
    - double mass    : mass, divided by the solar mass, assigned to _mass
    - char const* id : identifier to use when plotting the orbits
  */
  CelestialBody(vec3, vec3, double,char const*);

  /*
  double
  getMass()

  Returns the mass of the body
  Output:
    - _mass
  */
  double getMass();

  /*
  vec3
  &position()

  Output:
    - &_position
  */
  vec3 &position();

  /*
  vec3
  &velocity()

  Output:
    - &_velocity
  */
  vec3 &velocity();

  /*
  vec3
  &interaction()

  Output:
    - &_interaction
  */
  vec3 &interaction();

  char const* _id;

private:
  double _mass;
  vec3 _position, _velocity, _interaction;

};
#endif
