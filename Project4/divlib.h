#ifndef DIVLIB_H
#define DIVLIB_H

#include <iostream>

void clean_matrix(double **,int);
/*
  void
  clean_matrix(double ** A,int dim)

  Deallocates memory of A.

  Input:
    - double ** A   : quadratic matrix to be deallocated
    - int dim       : length of A[0]
*/

void print_matr(double ** ,int );
/*
  void
  print_matr(double ** A,int n)

  Writes to Unix shell/command line prompt the matrix A.
  Input:
    - double ** A   : quadratic matrix to be printed
    - int dim       : length of A[0]
*/

double ** eye(int );
/*
  double **
  eye(int n)

  Creates and return an identity matrix of dimension n*n

  Input:
    - int n  
  Output:
    An n*n identity matrix
*/
#endif
