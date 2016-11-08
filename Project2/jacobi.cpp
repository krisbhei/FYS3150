#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <ctime>
#include <fstream>

#include "problem.h" //iostream included
#include "divlib.h"

using namespace std;

ofstream ofile;

//Functions that solves the described problem in the report
void jacobi_solver(double ** ,double ** , int , double tol = 1E-8);
/*
  void
  jacobi_solver(double ** A,double ** R, int n, double tol)

  Finds the eigenvalues and vectors of A using Jacobi's rotation method.
  The eigenvalues will be the diagonal elements of A after the function call.
  The respective eigenvectors,i, will be stored in R at column i.

  Input:
    - double ** A   : quadratic,symmeric matrix to find the eigenvectors and eigenvalues
    - double ** R   : matrix where the eigenvectors are stored
    - int n         : the length of A[0]
    - double tol    : tolerance for when the function consider
                      a value to be zero
*/

double find_max(double **, int* , int);
/*
  double
  find_max(double ** A, int* ind, int n)

  Finds the offdiagonal element of A that has the largest
  absolute value.
  The function assigns ind the column and the row index
  of where the largest (in absolute value) element of A is located.

  Input:
    - double ** A   : quadratic,symmetric matrix to search for the offdiagonal element that has the largest
                      absolute value
    - double * ind  : array with length two (assumed) to store where the largest offdiagonal
                      element is located in A
    - int n         : the length of A[0]

  Output:
    The largest aboslute value of an offdiagonal element in A
*/

int * lowest_eig_ind(double ** ,int,int eig_n = 3);
/*
  int *
  lowest_eig_ind(double ** A,int dim)

  Assumed to be called after A has been undergone a similarity tranformation.
  Finds the indicies of where the three lowest eigenvalues are stored.
  The indicies are sorted such that the index at 0 gives the lowest eigenvalue,
  1 the second lowest and 2 the third lowest.

  Input:
    - double ** A   : diagonal,quadratic matrix
    - int n         : the length of A[0]
    - int eig_n     : number of eigenvalues
  Output:
    Array of length 3 containing the indicies to the three lowest eigenvalues
    in increasing order.
*/

void write_to_file_inter(double* ,double , double * , double ,int );
/*
  void
  write_to_file_inter(double* eigvector,double lambda, double * rho, double w_r,int dim)

  Called in interacting_case to write computed values of jacobi_solver.
  Here the probability density is computed and stored to file at the same row
  as the corresping value of rho_i.
  The file is written for readPlotFile.py to interpret and plot.

  Input:
    - double * eigvector  : resulting value of the wavefunction
    - double lambda       : eigenvalue to eigvector
    - double * rho
    - double w_r          : oscilator potential used to compute the
                            values of eigvector
    - int dim             : number of rows in eigvector
*/

//End: functions that solves the described problem
//------------------------------------------------

//Functions that are called by project2.py
extern "C" double test_find_max(double **,int);
/*
  double
  test_find_max(double ** A,int dim)

  Finding the max aboslute value of the offdiagonal
  element of A, and return the value of it.
  The wrapper file, project2.py compares the value
  with an expected value.

  Input:
    - double ** A   : quadratic,symmetric matrix to search for the offdiagonal element that has the largest
                      absolute value
    - int dim       : length of A[0]

  Output:
    The largest aboslute value of an offdiagonal element in A

*/

extern "C" int test_find_lowest(double**,int *,int);
/*
  int
  test_find_lowest(double** A,int * expected_ind,int dim)

  The function tests if the three indicies corresponding to the three
  lowest eigenvalues of A found in lowest_eig_ind does not differ
  greater than 1E-14 to the expected values given as arguments.
  The matrix A is assumed to be diagonal, wheras the diagonal
  elements are the 'eigenvalues'.
  project2.py interprets the int value returned by this function.
  If the calculation went wrong, the wrapperfile raises an AssertionError.

  Input:
    - double ** A         : diagonal,quadratic matrix
    - int * expected_ind  : array containg the expected indicies
    - int dim             : length of A[0]

  Output:
    An integer to state if the test went fine or not
    (0 for fine and -1 for incorrect indicies).
*/
extern "C" int test_solver(double **,double **,double*,int);
/*
  int
  test_solver(double ** A,double ** expected_eigvec,double * expected_eig,int dim)

  Tests if jacobi_solver gives expected eigenvalues and eigenvectors.
  project2.py interprets the int value returned by this function.
  The test_find_lowest is assumed to be called first.
  It is then possible to use lowest_eig_ind to sort and compare the expected eigenvalues in
  increasing order.

  If the calculation went wrong, the wrapperfile raises an AssertionError.

  Input:
    - double ** A               : quadratic matrix to find eigenvalues and eigenvectors for
    - double ** expected_eigvec : expected eigenvectors.
                                  column i corresponds to eigenvalue at index i in expected_eig.
    - double * expected_eig     : array containg the eigenvectors in increasing order.
    - int dim                   : dimension of A[0] (A is assumed quadratic)

  Output:
    An integer to state if the test went fine or not.
    -1 states that the computed eigenvectors are not normalized
    -2 states that the difference between the expected eigenvectors and the computed ones
       differ more than 1E-10
    -3 states that the difference between the expected eigenvalues and the computed ones
       differ more than 1E-10
*/

extern "C" double * non_interacting(double ,int);
/*
  double *
  non_interacting(double rho_max,int n)

  Sets up the system for the non-interacting case described in the report,
  and then solves it by calling jacobi_solver.

  Input:
    - double rho_max
    - int n           : number of grid points
  Output:
    The three lowest eigenvalues.
    They are used in project2.py to generate a tabular to
    show how the precision varies with n.
*/

extern "C" double non_interacting_time(double ,int);
/*
  double
  non_interacting_time(double rho_max,int n)

  Sets up the system for the non-interacting case described in the report,
  and then solves it by calling jacobi_solver. Here the execution time
  is also recorded to compare with numpy.linalg.eig in project2.py

  Input:
    - double rho_max
    - int n           : number of grid points
  Output:
    Total excecution time in seconds
*/

extern "C" void interacting_case(double, double, double *,int,int);
/*
 void
 interacting_case(double rho_max, double rho_max_1, double * w,int n,int w_n)

 Sets up the system for the interacting case described in the report,
 and then solves it by calling jacobi_solver.
 The result is written by calling write_to_file_inter to a file named: jacobi_lowest.dat

 The reason that the function takes in two values for rho_max is
 because of for small values (< 1) of the oscilator potential w
 cuts the wavefunction in the plots.

 Input:
   - double rho_max   : first max value of rho.
                        assumed to be used when the values of w are less than 1
   - double rho_max_1 : assumed to be used with values of w that are greater or equal to 1
   - double * w       : oscilator potentials
   - int n            : number of grid points
   - int w_n          : length of w
*/

extern "C" void interacting_case_noCoul(double, double, double *,int,int);
/*
 void
 interacting_case(double rho_max, double rho_max_1, double * w,int n,int w_n)

 Sets up the system for the interacting case described in the report,
 and then solves it by calling jacobi_solver.
 The result is written by calling write_to_file_inter to a file named: jacobi_lowest_noCoul.dat

 The reason that the function takes in two values for rho_max is
 because of for small values (< 1) of the oscilator potential w
 cuts the wavefunction in the plots.

 Input:
   - double rho_max   : first max value of rho.
                        assumed to be used when the values of w are less than 1
   - double rho_max_1 : assumed to be used with values of w that are greater or equal to 1
   - double * w       : oscilator potentials
   - int n            : number of grid points
   - int w_n          : length of w
*/

//End: functions that are called by project2.py
//---------------------------------------------

int main(int argc, char* argv []) //default call; solve for non interacting
{
  const double rho0 = 0;

  int n;
  double rho_max,rho_max_1;

  if (argc < 2)
  {
    n = 200;
    rho_max = 50;
    rho_max = 5;
  }
  else
  {
    n = atoi(argv[1]);
    rho_max = atof(argv[2]);
    rho_max_1 = atof(argv[3]);
  }

  int dim = n-2;
  const double u0 = 0, u_inf = 0;

  double ** eigvectors = eye(dim);
  double ** A = eye(dim);

  double * rho = generate_rho(rho0,rho_max,n);

  generate_system(A,rho,&V,n-2,0.01);
  jacobi_solver(A,eigvectors,n-2);
  return 0;
}

void jacobi_solver(double ** A,double ** R, int n, double tol)
{
  int iter = 0;
  int * ind = new int[2];

  while(find_max(A,ind,n) > tol)
  {
    int k = ind[0], l = ind[1];

    double t = 0, c = 0,s = 0;

    double tau = (A[l][l] - A[k][k])/(2.*A[k][l]);
    t = (tau > 0 ? 1./(tau + sqrt(1. + tau*tau)) : 1./(tau - sqrt(1. + tau*tau)));

    c = 1./sqrt(1+t*t);
    s = t*c;

    double a_kk = A[k][k];
    double a_ll = A[l][l];
    double a_kl = A[k][l];

    A[k][k] = a_kk*c*c - 2*a_kl*s*c + a_ll*s*s;
    A[l][l] = a_ll*c*c + 2*a_kl*s*c + a_kk*s*s;
    A[k][l] = 0.;
    A[l][k] = 0.;

    for (size_t i = 0; i < n; i++)
    {
      if (i != k && i != l)
      {
        double a_ik = A[i][k];
        double a_il = A[i][l];

        A[i][k] = a_ik*c - a_il*s;
        A[i][l] = a_il*c + a_ik*s;

        A[k][i] = A[i][k];
        A[l][i] = A[i][l];

      }
      double r_ik = R[i][k];
      double r_il = R[i][l];
      R[i][k] = c*r_ik - s*r_il;
      R[i][l] = c*r_il + s*r_ik;
    }
    iter++;
  }
  cout << "Number of iterations: "<< iter << endl;
  delete [] ind;
  return;
}//end: jacobi_solver

double find_max(double ** A, int* ind, int n)
{
  double max = 0;
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = i+1; j < n; j++)
    {
      if (fabs(A[i][j]) > max )
      {
        ind[0] = i;
        ind[1] = j;
        max = fabs(A[i][j]);
      }
    }
  }
 return max;
} //end: find_max

int * lowest_eig_ind(double ** A,int dim,int eig_n)
{
  int * ind = new int[eig_n];
  for (size_t i = 0; i < eig_n; i++)
  {
    ind[i] = i;
  }

  //start: sort
  int k = 1;
  while (k<eig_n)
  {
    int k1 = ind[k-1];
    int k2 = ind[k];

    double eig_prev = A[k1][k1];
    double eig = A[k2][k2];

    if (eig_prev > eig) //swap
    {
      ind[k-1] = k2;
      ind[k] = k1;
      k = 1;
    }
    else
    {
      k++;
    }
  }
  //end: sort

  for (size_t i = eig_n; i < dim; i++)
  {
    if (A[i][i] < A[ind[0]][ind[0]])
    {
      ind[2] = ind[1];
      ind[1] = ind[0];
      ind[0] = i;
    }
    else if (A[i][i] < A[ind[1]][ind[1]])
    {
      ind[2] = ind[1];
      ind[1] = i;
    }
    else if (A[i][i] < A[ind[2]][ind[2]])
    {
      ind[2] = i;
    }
  }
  return ind;
} //end: lowest_eig_ind

void write_to_file_inter(double* eigvector,double lambda, double * rho, double w_r,int dim)
{
  ofile << "omega_r = " << w_r << endl;
  ofile << "lambda = " << lambda << endl;
  for (size_t i = 0; i < dim; i++)
  {
    ofile << rho[i] << setw(20) << eigvector[i]*eigvector[i] << endl;
  }
  return;
} //end: write_to_file_inter



extern "C" double test_find_max(double ** A,int dim)
{
  int * ind = new int[2];
  return find_max(A,ind,dim);
} //end: test_find_max

extern "C" int test_find_lowest(double** A,int * expected_ind,int dim)
{
  int * ind = lowest_eig_ind(A,dim);
  for (size_t i = 0; i < 3; i++)
  {
    int i_ = ind[i];
    int ii_ = expected_ind[i];

    if (fabs(A[i_][i_] - A[ii_][ii_])>1E-14)
    {
      return -1;
    }
  }
  delete [] ind;
  return 0;
} //end: test_find_lowest

extern "C" int test_solver(double ** A,double ** expected_eigvec,double * expected_eig,int dim)
{
  double ** R = eye(dim);
  jacobi_solver(A,R,dim);
  int * ind = lowest_eig_ind(A,dim);

  const double tol = 1E-10;
  for (size_t i = 0; i < 3; i++)
  {
    int i_ = ind[i];
    if (fabs(A[i_][i_]  - expected_eig[i] )>tol)
    {
      return -3;
    }
    double sum = 0;
    for (size_t j = 0; j < dim; j++)
    {
      if (fabs(expected_eigvec[j][i] - R[j][i_]) > tol)
      {
        cout << expected_eigvec[j][i] << endl;
        cout <<  R[j][i_] << endl;
        return -2;
      }
      sum += R[j][i]*R[j][i];
    }

    if ((1.-sum) > tol)
    {
      return -1;
    }
  }

  clean_matrix(R,dim);
  return 0;
} //end:test_solver


extern "C" double * non_interacting(double rho_max,int n)
{
  const double rho0 = 0;

  int dim = n-2;


  double ** eigvectors = eye(dim);
  double ** A = eye(dim);

  double * rho = generate_rho(rho0,rho_max,n);

  generate_system(A,rho,&V,dim);

  jacobi_solver(A,eigvectors,dim);

  int * ind = lowest_eig_ind(A,dim);

  double * lowest_eig = new double[3];
  for (size_t i = 0; i < 3; i++)
  {
    lowest_eig[i] = A[ind[i]][ind[i]];
  }

  clean_matrix(A,dim);
  clean_matrix(eigvectors,dim);
  delete [] rho;

  cout << "Done solving for non-interacting case with n = "<<n;
  cout << " and rho_max = " <<rho_max <<endl;
  return lowest_eig;
} //end: non_interacting

extern "C" double non_interacting_time(double rho_max,int n)
{
  const double rho0 = 0;

  //n-2 since we skip the boundaries
  int dim = n-2;

  double ** eigvectors = eye(dim);
  double ** A = eye(dim);

  double * rho = generate_rho(rho0,rho_max,n);

  generate_system(A,rho,&V,dim);

  double start = clock();
  jacobi_solver(A,eigvectors,dim);
  double end = clock();

  clean_matrix(A,dim);
  clean_matrix(eigvectors,dim);
  delete [] rho;

  return (end-start)/CLOCKS_PER_SEC;
} //end: non_interacting_time

extern "C" void interacting_case(double rho_max, double rho_max_1, double * w,int n,int w_n)
{
  const double rho0 = 0;
  const double u0 = 0, u_inf = 0;
  int dim = n-2;

  double ** eigvectors = eye(dim);
  double ** A = eye(dim);

  double * rho = generate_rho(rho0,rho_max,n);

  ofile.open("jacobi_lowest.dat");
  ofile << "N = " << n << endl;
  for (size_t i = 0; i < w_n; i++)
  {
    if (w[i] >= 1.)
    {
      rho = generate_rho(0,rho_max_1,n);
    }
    generate_system(A,rho,&V_inter,dim,w[i]);

    jacobi_solver(A,eigvectors,dim);

    int ind = lowest_eig_ind(A,dim)[0];

    double * eigvector = new double[n];
    eigvector[0] = u0;
    for (size_t i = 0; i < dim; i++)
    {
      eigvector[i+1] = eigvectors[i][ind];
    }
    eigvector[n-1] = u_inf;

    write_to_file_inter(eigvector,A[ind][ind],rho,w[i],n);

    eigvectors = eye(dim);
    cout << "Done solving for interacting case with omega_r = " << w[i] << endl;
  }
  ofile.close();

  clean_matrix(A,dim);
  clean_matrix(eigvectors,dim);
  delete [] rho;

  return;
}//end: interacting_case

extern "C" void interacting_case_noCoul(double rho_max, double rho_max_1, double * w,int n,int w_n)
{
  const double rho0 = 0;
  const double u0 = 0, u_inf = 0;
  int dim = n-2;

  double ** eigvectors = eye(dim);
  double ** A = eye(dim);

  double * rho = generate_rho(rho0,rho_max,n);

  ofile.open("jacobi_lowest_noCoul.dat");
  ofile << "N = " << n << endl;
  for (size_t i = 0; i < w_n; i++)
  {
    if (w[i] >= 1.)
    {
      rho = generate_rho(0,rho_max_1,n);
    }
    generate_system(A,rho,&V_inter_nCoul,dim,w[i]);

    jacobi_solver(A,eigvectors,dim);

    int ind = lowest_eig_ind(A,dim)[0];

    double * eigvector = new double[n];
    eigvector[0] = u0;
    for (size_t i = 0; i < dim; i++)
    {
      eigvector[i+1] = eigvectors[i][ind];
    }
    eigvector[n-1] = u_inf;

    write_to_file_inter(eigvector,A[ind][ind],rho,w[i],n);

    eigvectors = eye(dim);
    cout << "Done solving for interacting case with omega_r = " << w[i] << endl;
  }
  ofile.close();

  clean_matrix(A,dim);
  clean_matrix(eigvectors,dim);
  delete [] rho;

  return;
}//end: interacting_case
