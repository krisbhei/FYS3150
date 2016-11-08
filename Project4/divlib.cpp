#include "divlib.h"

int periodic(int N,int index)
{
  return index == N ? 0 : index; //0 til N-1
}
void clean_matrix(double ** A,int dim)
{
  for (int i = 0; i < dim; i++)
  {
    delete [] A[i];
  }
  delete [] A;
}
void print_matr(double ** A,int n)
{
  std::cout << "---------------" << std::endl;
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
    printf("%7.2f",A[i][j]);
    }
    printf("\n");
  }
  std::cout << "---------------" << std::endl;
}

double ** eye(int n)
{
  double ** mat = new double*[n];

  for (size_t i = 0; i < n; i++)
  {
    mat[i] = new double[n];
    for (size_t j = 0; j < n; j++)
    {
      mat[i][j] = 0;
    }
    mat[i][i] = 1;
  }
  return mat;
}
