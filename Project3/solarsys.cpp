#include "solarsys.h"
#include "differential.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>

#define G 4*M_PI*M_PI
using namespace std;

SolarSystem::SolarSystem()
{
  _celestialBodies.push_back(SUN);
  _numBodies = 1;
  _kineticEnergy = 0.;
  _potentialEnergy = 0.;
  _momentum.zeros();
}

void SolarSystem::addCelestialBody(CelestialBody body)
{
  _celestialBodies.push_back(body);
  ++_numBodies;
}

void SolarSystem::setCenterOfMassOrigin()
{
  double totalMass = 0;
  _centerOfMass.zeros();

  for (CelestialBody &body: _celestialBodies)
  {
    totalMass += body.getMass();
    _centerOfMass += body.position()*body.getMass();
  }

  _centerOfMass /= totalMass;
  for (CelestialBody &body : _celestialBodies)
  {
    body.position() -= _centerOfMass;
  }
}

void SolarSystem::setInitialVelocitySun()
{
  vec3 zeroMomentumVelocity;

  for(CelestialBody &body : _celestialBodies)
  {
    zeroMomentumVelocity -= body.getMass()*body.velocity();
  }
  _celestialBodies[0].velocity() = zeroMomentumVelocity;
  _celestialBodies[0]._id = "Sun";
}

std::vector<CelestialBody> &SolarSystem :: celestialBodies()
{
  return _celestialBodies;
}

vec3 SolarSystem::newtonianGravitationalForce(CelestialBody body1, CelestialBody body2)
{
  vec3 pos1 = body1.position();
  vec3 pos2 = body2.position();

  double dr = (pos1-pos2).length();
  double dr3 = dr*dr*dr;
  double gm = -G/dr3;

  vec3 acc = gm*(pos1-pos2);
  return acc;
}
vec3 SolarSystem::relativisticCorrection(CelestialBody body1, CelestialBody body2)
{
  vec3 acc = newtonianGravitationalForce(body1,body2);

  vec3 pos_body1 = body1.position();
  vec3 pos_body2 = body2.position();
  vec3 r = pos_body1 - pos_body2;

  double dist_squared = r.lengthSquared();
  double l_squared = (r.cross(body1.velocity())).lengthSquared();
  double c = 63198; //AU per year
  return acc*(1. + (3*l_squared)/(dist_squared*c*c));
}

void SolarSystem::computeInteractionRelativistic()
{
  for(CelestialBody &body : _celestialBodies)
  {
    body.interaction().zeros();
  }
  _momentum.zeros();
  _angularMomentum.zeros();
  _kineticEnergy = 0;
  _potentialEnergy = 0;
  for (size_t i = 0; i < _numBodies; i++)
  {
    CelestialBody &body1 = _celestialBodies[i];

    for (size_t j = i+1; j < _numBodies; j++)
    {
      CelestialBody &body2 = _celestialBodies[j];
      vec3 correction = relativisticCorrection(body1,body2);

      body1.interaction() += correction*body2.getMass();
      body2.interaction() -= correction*body1.getMass();
      _potentialEnergy += -(G*body2.getMass()*body1.getMass())/(body1.position() - body2.position()).length();

    }

    _kineticEnergy += .5*body1.getMass()*body1.velocity().lengthSquared();
  }
  return;
}



void SolarSystem::computeInteraction()
{
  for(CelestialBody &body : _celestialBodies)
  {
    body.interaction().zeros();
  }
  _momentum.zeros();
  _angularMomentum.zeros();
  _kineticEnergy = 0;
  _potentialEnergy = 0;

  for (size_t i = 0; i < _numBodies; i++)
  {
    CelestialBody &body1 = _celestialBodies[i];

    for (size_t j = i+1; j < _numBodies; j++)
    {
      CelestialBody &body2 = _celestialBodies[j];

      vec3 force = newtonianGravitationalForce(body1,body2);

      body1.interaction() += force*body2.getMass();
      body2.interaction() -= force*body1.getMass();
      _potentialEnergy += -(G*body2.getMass()*body1.getMass())/(body1.position() - body2.position()).length();
    }

    double mass_body1 = body1.getMass();
    vec3 velocity_body1 = body1.velocity();
    vec3 position_body1 = body1.position();

    _momentum += mass_body1 * velocity_body1;
    _kineticEnergy += .5*mass_body1 * velocity_body1.lengthSquared();
    _angularMomentum += position_body1.cross( (velocity_body1 * mass_body1) );
  }
  return;
}

template<typename T>
void SolarSystem::takeTimeSolver(T odeSolver,double t_min,double t_max,int n,double &totaltime)
{
  odeSolver.setTimeStep((t_max-t_min)/(n-1.));
  double start = clock();
  for (int i = 0; i < n; i++)
  {
    odeSolver.calculateOneStep(this);
  }
  double end = clock();
  totaltime = double(end-start)/CLOCKS_PER_SEC;
  return;
}
void SolarSystem::computeMovementRelativistic(double t_min,double t_max,long n,string filename)
{
  double thetaPrev = 0;
  double theta = 0;

  double rPreviousPrevious 	= 0;
  double rPrevious   	 	= 0;
  double r 		 	= 0;

  VerletRelativistic odeSolver;
  double timeStep = (t_max - t_min)/(n-1.);
  odeSolver.setTimeStep(timeStep);

  int j = 0;

  vec3 previousPosition(0,0,0);

  _file.open(filename,ofstream::out);
  for (size_t i = 0; i < n; i++)
  {
    odeSolver.calculateOneStep(this);
    CelestialBody mercury = _celestialBodies[1];
    CelestialBody sun = _celestialBodies[0];

    double rCurrent = (mercury.position()-sun.position()).length();
    if ( rCurrent > rPrevious && rPrevious < rPreviousPrevious )
    {
        double x = previousPosition.x();
		    double y = previousPosition.y();

        _file << "Time: " << (t_min + double(i)*timeStep) << endl;
        _file << "Perihelion: " << atan2(y,x)*(3600*180)/M_PI << endl;
        _file << "Total energy: " <<  _kineticEnergy + _potentialEnergy << endl;
        ++j;
	  }
    rPreviousPrevious 	= rPrevious;
	  rPrevious		= rCurrent;

    previousPosition	= mercury.position() - sun.position();
  }
  _file.close();
}

template<typename T>
void SolarSystem::computeMovement(T odeSolver,double t_min,double t_max,int n,string filename)
{

  _file.open(filename,ofstream::out);
  _file << "Center of mass: " << _centerOfMass.x() << " " << _centerOfMass.y() << " "<< _centerOfMass.z() << endl;
  _file << "Number of bodies: " <<_numBodies<<endl;
  _file << "Gridpoints: " << n << endl;
  _file <<"Identifier(s): ";
  for (size_t i = 0; i < _numBodies; i++)
  {
    _file << _celestialBodies[i]._id<< " ";
  }
  _file << endl;

  double h = (t_max-t_min)/(n-1.);
  double current_year = t_min;
  odeSolver.setTimeStep(h);

  for (int i = 0; i < n; i++)
  {
    current_year = t_min + i*h;
    _file << "Current year: "<< current_year << endl;

    odeSolver.calculateOneStep(this);

    for (int j = 0 ; j< _numBodies;j++)
    {
      vec3 pos = _celestialBodies[j].position();
      vec3 vel = _celestialBodies[j].velocity();

      _file <<"location: "<< setprecision(15)<< pos.x() << " " << pos.y()<< " "<<pos.z() <<endl;
      _file <<"velocity: "<< setprecision(15)<< vel.x() << " " << vel.y()<< " "<<vel.z() <<endl;
    }

    _file << "-----" << endl;
    _file << "total kinetic energy of system: " << _kineticEnergy << endl;
    _file << "total potential energy of system: " << _potentialEnergy << endl;
    _file << "angular momentum: " << _angularMomentum.x();
    _file << " " << _angularMomentum.y();
    _file << " " << _angularMomentum.z() << endl;
    _file << "momentum: " << _momentum.x();
    _file << " " << _momentum.y();
    _file << " " << _momentum.z() << endl;
    _file << "-----" << endl;

  }
  _file.close();
  return;
}
void SolarSystem::reset()
{
  _celestialBodies.clear();
  SUN.velocity().zeros();
  SUN.position().zeros();
  SUN.interaction().zeros();
  SUN._id = "-";
  _celestialBodies.push_back(SUN);
  _numBodies = 1;
}


template void SolarSystem::takeTimeSolver<Euler>(Euler odeSolver,double t_min,double t_max,int n,double &totaltime);
template void SolarSystem::takeTimeSolver<Verlet>(Verlet odeSolver,double t_min,double t_max,int n,double &totaltime);

template void SolarSystem::computeMovement<Euler>(Euler odeSolver,double t_min,double t_max,int n,string filename);
template void SolarSystem::computeMovement<Verlet>(Verlet odeSolver,double t_min,double t_max,int n,string filename);
