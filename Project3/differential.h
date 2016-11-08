#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H
#include "vec3.h"
#include "celestialbody.h"
#include "solarsys.h"


class Euler
{
public:
  void calculateOneStep(SolarSystem*);
  void setTimeStep(double);
private:
  double _dt;
};

class VerletRelativistic
{
public:
  void calculateOneStep(SolarSystem*);
  void setTimeStep(double);
private:
  double _dt;
};

class Verlet
{
public:
  void calculateOneStep(SolarSystem*);
  void setTimeStep(double);
private:
  double _dt;
};

#endif
