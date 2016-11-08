#ifndef PROBLEM_H
#define PROBLEM_H

void generate_system(double**,double * , double (*potential)(double,double),\
                     int,double w_r = 0);
/*
  void
  generate_system(double**A,double * rho, double (*potential)(double,double),
                  int n,double w_r)

  Generates the systems to solve as described in the report.

  Input:
    - double ** A                        : matrix to store the values of the system.
    - double * rho                       : array of length n containing the values
                                           of rho in the report
    - double (*potential)(double,double) : function describing the potential
                                           to be used in the system.
    - int n                              : number of grid points
    - double w_r                         : oscilator potential (if there is any)
*/
double * generate_rho(double, double , int );
/*
  double *
  generate_rho(double rho0, double rhoN, int dim)

  Generates dim uniformly spaced values within the intervall [rho0,rhoN].

  Input:
    - double rho0
    - double rhoN
    - int n         : number of grid points
  Output:
    Array containing the n uniformly spaced values of rho
*/

double V(double,double w_r=1);
/*
  double
  V(double rho,double w_r)

  The potential used in the non-interacting case.
  Input:
    - double rho
    - double w_r  : dummy argument to define the potential functions as arguments
                    to generate_system
  Output:
    rho*rho
*/

double V_inter(double , double );
/*
  double
  V(double rho,double w_r)

  The potential used in the interacting case.
  Input:
    - double rho
    - double w_r  : oscilator potential
  Output:
    rho*rho*w_r*w_r + 1./rho
*/
double V_inter_nCoul(double, double);
/*
  double
  V(double rho,double w_r)

  The potential used in the interacting case without Coulomb
  interaction.
  Input:
    - double rho
    - double w_r  : oscilator potential
  Output:
    rho*rho*w_r*w_r + 1./rho
*/
#endif
