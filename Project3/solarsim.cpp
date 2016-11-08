/*
These functions are called by the wrapper file project3.py
*/

#include "celestialbody.h"
#include "solarsys.h"
#include "differential.h"
#include <iostream>
#include <cmath>
#include <string>

using namespace std;
/*
Computes the perihelion using relativistic correction
*/
extern "C" void generalRelativity(double min_yr,double max_yr,const char* filename, int n)
{
  SolarSystem solarSystem;
  CelestialBody mercury(vec3(.3075,0,0),vec3(0,12.44,0),1.2*1E-7,"Mercury");
  solarSystem.addCelestialBody(mercury);
  solarSystem.setCenterOfMassOrigin();
  solarSystem.computeMovementRelativistic(min_yr,max_yr,n,string(filename));
  return; //End: generalRelativity
}

/*
Runs the simulation for Sun-Earth-Jupiter to check for zero momentum, then runs for all planets
with initial positions and velocities from NASA
*/
extern "C" void allPlanets(double min_yr,double max_yr,double max_yr_three, const char* filename, int n_all,int n_three)
{
  SolarSystem solarSystem;
  string file = string(filename);

  //Start: Run for the three-body system
  CelestialBody earth(vec3(1,0,0),vec3(0,2*M_PI,0),3*1E-6,"Earth");
  CelestialBody jupiter(vec3(5.2,0,0),vec3(0,M_PI,0),.95*1E-3,"Jupiter");

  solarSystem.addCelestialBody(jupiter);
  solarSystem.addCelestialBody(earth);

  solarSystem.setCenterOfMassOrigin();
  solarSystem.setInitialVelocitySun();

  string f = string("threeBodyZeroMomentum_")+file;
  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr_three,n_three,f);
  //End: Run for the three-body system

  //Start: Run for all planets in the solar system (including Pluto) with data from NASA
  //Date for the planet's initial position and velocity: 2016-10-05
  solarSystem.reset();

  CelestialBody NASA_earth(vec3(9.819128739328793E-01,2.104822076393571E-01,0),vec3(-3.851159854117840E-03*365,1.677807321756382E-02*365,0),3*1E-6,"Earth");
  CelestialBody NASA_mercury(vec3(-1.388351215994794E-01,2.874076124640064E-01,0),vec3(-3.081033504804020E-02*365,-1.153752302730325E-02 *365,0),1.2*1E-7,"Mercury");
  CelestialBody NASA_venus(vec3(2.531171709934956E-03,-7.235439738148072E-01,0),vec3(2.008934322426700E-02*365,-9.726053625259401E-05*365,0),2.45*1E-6,"Venus");
  CelestialBody NASA_mars(vec3(1.074137801923908E+00,-8.751565791508180E-01,0),vec3(9.406114898320066E-03*365,1.202410268499505E-02*365,0),3.3*1E-7,"Mars");
  CelestialBody NASA_jupiter(vec3(-5.433468170028908E+00,-3.819061221110369E-01,0),vec3(4.425651679847022E-04*365,-7.171108917491057E-03*365,0),.95*1E-3,"Jupiter");
  CelestialBody NASA_saturn(vec3(-2.318303159285024E+00,-9.761896118742531E+00,0),vec3(5.122777078109817E-03*365,-1.306326757884738E-03*365,0),2.75*1E-4,"Saturn");
  CelestialBody NASA_uranus(vec3(1.847838443295973E+01,7.526847462019028E+00,0),vec3(-1.512383759608680E-03*365,3.459146288519939E-03*365,0),4.4*1E-5,"Uranus");
  CelestialBody NASA_neptune(vec3(2.825072704568992E+01, -9.952093577677799E+00,0),vec3(1.022588623946866E-03*365, 2.979574756810357E-03*365,0),.515*1E-4,"Neptune");
  CelestialBody NASA_pluto(vec3(9.393096450667111E+00, -3.182064102580347E+01,0),vec3(3.065499934972441E-03*365, 2.293283900283695E-04*365,0),.655*1E-8,"Pluto");

  solarSystem.addCelestialBody(NASA_earth);
  solarSystem.addCelestialBody(NASA_mercury);
  solarSystem.addCelestialBody(NASA_venus);
  solarSystem.addCelestialBody(NASA_mars);
  solarSystem.addCelestialBody(NASA_jupiter);
  solarSystem.addCelestialBody(NASA_saturn);
  solarSystem.addCelestialBody(NASA_uranus);
  solarSystem.addCelestialBody(NASA_neptune);
  solarSystem.addCelestialBody(NASA_pluto);

  solarSystem.setCenterOfMassOrigin();
  solarSystem.setInitialVelocitySun();

  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr,n_all,file);
  //End: Run for all planets in the solar system (including Pluto) with data from NASA

  return; //End: allPlanets
}
extern "C" void sunEarthJupiter(double inc_factor,double min_yr,double max_yr,const char* filename, int n)
{
  SolarSystem solarSystem;
  CelestialBody earth(vec3(1,0,0),vec3(0,2*M_PI,0),3*1E-6,"Earth");
  CelestialBody jupiter(vec3(5.2,0,0),vec3(0,M_PI,0),.95*1E-3*inc_factor,"Jupiter");

  solarSystem.addCelestialBody(jupiter);
  solarSystem.addCelestialBody(earth);

  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr,n,string(filename));
  return;
}
extern "C" void escapeSun(double vx,double vy,double min_yr,double max_yr,const char* filename, int n)
{
  SolarSystem solarSystem;
  CelestialBody planet(vec3(1,0,0),vec3(vx,vy,0),3*1E-6,"Planet");

  //Start: Compute orbit for the experimental velocity
  solarSystem.addCelestialBody(planet);

  string file = string(filename);
  string file1 = string("velocity_experiment_") + file;
  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr,n,file1);
  //End: Compute orbit for the experimental velocity

  //Start: Compute orbit with excat escaping velocity
  double exact_velocity_escape = sqrt(8*M_PI*M_PI);
  string file2 = string("velocity_exact_") + file;

  solarSystem.reset();

  CelestialBody planet_ex(vec3(1,0,0),vec3(0,exact_velocity_escape,0),3*1E-6,"Planet");
  solarSystem.addCelestialBody(planet_ex);

  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr,n,file2);
  //End: Compute orbit with excat escaping velocity

  return; //End: escapeSun
}
extern "C" double timeVerlet(double min_yr,double max_yr, int n)
{
  SolarSystem solarSystem;

  CelestialBody earth(vec3(1,0,0),vec3(0,2*M_PI,0),3*1E-6,"Earth");
  solarSystem.addCelestialBody(earth);

  double finalTime = 0;
  solarSystem.takeTimeSolver(Verlet(),min_yr,max_yr,n,finalTime);
  return finalTime;
}
extern "C" double timeEuler(double min_yr,double max_yr, int n)
{
  SolarSystem solarSystem;

  CelestialBody earth(vec3(1,0,0),vec3(0,2*M_PI,0),3*1E-6,"Earth");
  solarSystem.addCelestialBody(earth);

  double finalTime = 0;
  solarSystem.takeTimeSolver(Euler(),min_yr,max_yr,n,finalTime);
  return finalTime;
}
extern "C" void testEarthSun(double min_yr,double max_yr,const char* filename, int n)
{
  SolarSystem solarSystem;
  CelestialBody earth(vec3(1,0,0),vec3(0,2*M_PI,0),3*1E-6,"Earth");
  solarSystem.addCelestialBody(earth);

  string file_e = string("e_")+string(filename);
  string file_v = string("v_")+string(filename);

  solarSystem.computeMovement<Verlet>(Verlet(),min_yr,max_yr,n,file_v);
  solarSystem.computeMovement<Euler>(Euler(),min_yr,max_yr,n,file_e);
  return; //End: testEarthSun
}

int main()
{
  testEarthSun(0,1,"earth-sun.txt",1E4);
  return 0;
}
