#ifndef JACOBI_H
#define JACOBI_H

double find_max(double **, int* , int);
void jacobi_solver(double ** ,double ** , double , double tol = 1E-8);
int * lowest_eig_ind(double ** ,int );
double test_find_max(double ** );
void write_to_file_inter(double* ,double , double * , double ,int );

extern "C" double * non_interacting(double ,int);
extern "C" double non_interacting_time(double ,int);
extern "C" void interacting_case(double, double, double *,int,int);
#endif
