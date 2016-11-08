#include "problem.h"

void generate_system(double**A,double * rho, double (*potential)(double,double),int n,double w_r)
{
  double h = rho[1]-rho[0];

  for (int i = 1; i < n; i++)
  {
    A[i-1][i-1] = 2./(h*h) + potential(rho[i],w_r);
    A[i-1][i] = -1./(h*h);
    A[i][i-1] = A[i-1][i];
  }
  A[n-1][n-1] = 2./(h*h) + potential(rho[n],w_r);
}

double * generate_rho(double rho0, double rhoN, int dim)
{
    double h = (rhoN - rho0)/(dim-1);
    double * rho = new double[dim];

    for (int i = 0; i < dim; i++)
    {
      rho[i] = rho0 + i*h;
    }
    return rho;
}

double V(double rho,double w_r)
{
  return rho*rho;
}

double V_inter(double rho, double w_r)
{
  return (w_r*w_r)*(rho*rho) + 1./rho;
}
double V_inter_nCoul(double rho, double w_r)
{
  return (w_r*w_r)*(rho*rho);
}
