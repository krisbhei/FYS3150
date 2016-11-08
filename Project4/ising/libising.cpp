#include <stdlib.h>
#include <random>
#include "libising.h"

int periodic(int index, int len, int add)
{
    return (index+len+add)%len;

} //End: periodic

void initialize(int dim, int** spins,double & energy, double & magnetization)
{
    for (int i = 0; i<dim ; i++)
    {
        for (int j = 0; j<dim ; j++)
        {
            spins[i][j] = 1;
            magnetization += (double) spins[i][j];
        }
    }
    for (int i = 0; i<dim ; i++)
    {
        for (int j = 0; j<dim ; j++)
        {
            energy -= (double) spins[i][j]*(spins[periodic(i,dim,-1)][j] + spins[i][periodic(j,dim,-1)]);
        }
    }
    return;

} //End: initialize

void metropolis(int** spins,int dim, long trials, double T,double expectations[5],double energy, double magnetization)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distr(0.0,1.0);

    double w[17];
    for(int i = 0; i < 17 ; i++) w[i] = 0;
    for(int i = -8; i < 9 ; i+=4) w[i+8] = exp(-((double)i)/T);

    for(int cycle = 1 ; cycle <= trials ; cycle ++ )
    {
        for(int i = 0 ; i < dim*dim ; i++)
        {
            int r_x = (int) (distr(gen)*(double)dim);
            int r_y = (int) (distr(gen)*(double)dim);

            int deltaEnergy = 2*spins[r_y][r_x]*
                    (spins[r_y][periodic(r_x,dim,-1)] +
                    spins[periodic(r_y,dim,-1)][r_x] +
                    spins[r_y][periodic(r_x,dim,1)] +
                    spins[periodic(r_y,dim,1)][r_x]);

            if( distr(gen) <= w[deltaEnergy + 8])
            {
                spins[r_y][r_x] *= -1;
                magnetization += (double) 2*spins[r_y][r_x];
                energy += (double) deltaEnergy;
            }
        }

        expectations[0] += energy;
        expectations[1] += energy*energy;
        expectations[2] += magnetization;
        expectations[3] += magnetization*magnetization;
        expectations[4] += fabs(magnetization);
    }
    return;

} //End: metropolis

int ** init_matr(int dim)
{
    int ** matr = new int*[dim];
    for (int i = 0; i<dim; i++)
    {
        matr[i] = new int[dim];
        for (int j = 0; j<dim;j++)
        {
            matr[i][j] = 0;
        }
    }
    return matr;

} //End: init_matr

