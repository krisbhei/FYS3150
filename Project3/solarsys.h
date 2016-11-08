#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include <vector>
#include "celestialbody.h"
#include <fstream>
#include "vec3.h"
#include <string>

class SolarSystem
{
public:
  /*
  SolarSystem()

  Initializes the system. An CelestialBody which is supposed to be the Sun, that is
  SUN with position at zero and zero velocity is added at the end of _celestialBodies.
  Setting _numBodies to be 1, and sets _kineticEnergy, _potentialEnergy and _momentum
  to be zero.
  */
  SolarSystem();

  /*
  void
  addCelestialBody(CelestialBody body)

  Adds an object of type Celestial body into the system.
  The body is placed at the end in vector _celestialBodies.
  _numBodies is also incremented.

  Input:
    - CelestialBody body: The body to be stored
  */
  void addCelestialBody(CelestialBody);

  /*
  void
  setCenterOfMassOrigin()

  Sets the center of mass to be the origin.
  */
  void setCenterOfMassOrigin();

  /*
  void
  setInitialVelocitySun()

  Sets the velocity of the Sun such that the total momentum of the
  system is zero.
  Sets also an identifier to plot for the Sun's orbit.
  */
  void setInitialVelocitySun();

  /*
  void
  computeInteractionRelativistic()

  Computes the force between each body using the relativistic correction
  to the newtonian force
  */
  void computeInteractionRelativistic();

  /*
  void
  computeInteraction()

  Computes the force between each body using the newtonian force
  */
  void computeInteraction();

  /*
  vec3
  newtonianGravitationalForce(CelestialBody body1, CelestialBody body2)

  Computes a vector containing the force applied to body1 and body2.
  Input:
    - CelestialBody body1 : the body which has has the force vector pointing towards body2
    - CelestialBody body2 : body which has the force vector pointing at body1

  Output:
    - vec3 acc : the newtontial gravitional force to be applied to body1 and body2
  */
  vec3 newtonianGravitationalForce(CelestialBody, CelestialBody);

  /*
  vec3
  relativisticCorrection(CelestialBody body1, CelestialBody body2)

  Computes a the newtonian force vector by calling newtonianGravitationalForce and
  scales it with relativistic correction.
  Input:
    - CelestialBody body1 : the body which has has the force vector pointing towards body2
    - CelestialBody body2 : body which has the force vector pointing at body1

  Output:
    - vec3 : the newtontial gravitional force with relativistic correction to be applied to body1 and body2
  */
  vec3 relativisticCorrection(CelestialBody, CelestialBody);

  /*
  void
  computeMovement(T odeSolver,double t_min,double t_max,int n,string filename)

  Uses the given odeSolver to solve for the position and velocity for every body within the system.
  For every timestep, the position and velocity of every body is written to file.

  Input:
    - T odeSolver    : The solver class which conatins a specific intergration method
    - double t_min   : Start year for the bodies' orbit
    - double t_max   : End year for the bodies' orbit
    - int n          : Meshpoints
    - string filename:
  */
  template<typename T>
  void computeMovement(T odeSolver,double,double,int,std::string);

  /*
  void
  takeTimeSolver(T odeSolver,double t_min,double t_max,int n,string filename, double &totaltime)

  Uses the given odeSolver to solve for the position and velocity for every body within the system.
  Takes time for how long it takes for the program to solve the bodies' orbit.
  Sets totaltime to be the time taken.

  Input:
    - T odeSolver      : The solver class which conatins a specific intergration method
    - double t_min     : Start year for the bodies' orbit
    - double t_max     : End year for the bodies' orbit
    - int n            : Meshpoints
    - string filename  :
    - double &totaltime:
  */
  template<typename T>
  void takeTimeSolver(T odeSolver,double,double,int,double&);

  /* VENT
  void
  computeMovementRelativistic(double t_min,double t_max,int n,string filename)

  Uses an object of VerletRelativistic to solve for the position and velocity for every body within the system.
  Finds the perihelion angle and writes it to file.

  Input:
    - double t_min     : Start year for the bodies' orbit
    - double t_max     : End year for the bodies' orbit
    - long n            : Meshpoints
    - string filename  :
  */
  void computeMovementRelativistic(double,double,long,std::string);

  /*
  void
  reset()

  Clears all of the bodies in the system.
  */
  void reset();

  std::vector<CelestialBody> _celestialBodies;
  std::vector<CelestialBody> &celestialBodies();

private:
  CelestialBody SUN;
  int _numBodies;

  double _kineticEnergy;
  double _potentialEnergy;

  vec3 _centerOfMass;
  vec3 _momentum;
  vec3 _angularMomentum;

  std::ofstream _file;
};


#endif
