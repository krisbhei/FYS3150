#ifndef LIBISING_H
#define LIBISING_H

/* void
 * metropolis(int** spins,int dim, long trials, double T,double expectations[5],double energy, double magnetization)
 *
 * Runs the metropolis algorithm for a given number of trials. Computes the expectationvalues after each sweep thorugh the lattice.
 *
 * Input:
 *  - int ** spins           : configuration of the system
 *  - int dim                : number of spins
 *  - long trials            : number of total trials
 *  - double T               : temperature to calculate the energy
 *  - double expectations[5] : array to store the expectation values
 *  - double energy          : energy from the initial configuration
 *  - double magnetization   : magnetic moment from the inital configuration
*/
void metropolis(int**,int, long, double ,double expectations[5],double, double);

/*  int
 *  periodic(int index, int len, int add)
 *
 * Computes the index of an element in an array either in front(add=1) or the back(add=-1) of the element at index index.
 * The array is assumed to be len-periodic
 *
 * Input:
 *  - int index
 *  - int len
 *  - int add : either +1 or -1. Indicates if the program shall find the elememnt behind(-1) or in front (+1) of the element at index index
*/
int periodic(int, int, int);

/* int **
 * init_matr(int dim)
 *
 * Initializes a dim by dim matrix with elements 0.
 *
 * Input:
 *  - int dim
 *
 * Output:
 *  int ** matr : the initlized matrix
*/
int ** init_matr(int);

/* void
 * initialize(int dim, int** spins,double & energy, double & magnetization)
 *
 * Sets up the system being at the ground state with all spins pointing up
 *
 * Input:
 *  - int dim                : number of spins
 *  - int ** spins           : the lattice of spins
 *  - double & energy        : to store the energy of the initialized system
 *  - double & magnetization : to store the magnetic moment of the initialized system
*/
void initialize(int, int**,double &, double &);

/* void
 * initializeRandom(int dim, int** spins,double & energy, double & magnetization)
 *
 * Sets up the system being at a random state. The function calls the Mersenne Twister
 * for each spin and checks if the number is less or equal to .5.
 * If the check is true, then the current spin is set to point up (1).
 * If false, then the spin points down (-1).
 *
 * Input:
 *  - int dim                : number of spins
 *  - int ** spins           : the lattice of spins
 *  - double & energy        : to store the energy of the initialized system
 *  - double & magnetization : to store the magnetic moment of the initialized system
*/
void initializeRandom(int, int**,double &, double &);

#endif // LIBISING_H
