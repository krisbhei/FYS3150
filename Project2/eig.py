#This program is taken from the lecture notes with some modifications.
#See references in the report.
#Called by project2.py and compares time solving the system for the non-interacting case

#Program which solves the one-particle Schrodinger equation
#for a potential specified in function
#potential(). This example is for the harmonic oscillator in 3d


import numpy as np
import time

#Function for initialization of parameters
def initialize():
    RMin = 0.0
    RMax = 10.0
    lOrbital = 0
    Dim = 400
    return RMin, RMax, lOrbital, Dim

# Here we set up the harmonic oscillator potential
def potential(r):
    return r*r

def find_eig(RMin,RMax,lOrbital,Dim):
	#Initialize constants
	Step    = float(RMax)/(Dim+1)
	DiagConst = 2.0 / (Step*Step)
	NondiagConst =  -1.0 / (Step*Step)
	OrbitalFactor = lOrbital * (lOrbital + 1.0)

	#Calculate array of potential values
	v = np.zeros(Dim)
	r = np.zeros(Dim)
	for i in xrange(Dim):
		r[i] = RMin + (i+1) * Step;
		v[i] = potential(r[i]) + OrbitalFactor/(r[i]*r[i]);

	#Setting up a tridiagonal matrix and finding eigenvectors and eigenvalues
	Hamiltonian = np.zeros((Dim,Dim))
	Hamiltonian[0,0] = DiagConst + v[0];
	Hamiltonian[0,1] = NondiagConst;
	for i in xrange(1,Dim-1):
		Hamiltonian[i,i-1]  = NondiagConst;
		Hamiltonian[i,i]    = DiagConst + v[i];
		Hamiltonian[i,i+1]  = NondiagConst;
	Hamiltonian[Dim-1,Dim-2] = NondiagConst;
	Hamiltonian[Dim-1,Dim-1] = DiagConst + v[Dim-1];
	# diagonalize and obtain eigenvalues, not necessarily sorted

	starttime = time.time()
	EigValues, EigVectors = np.linalg.eig(Hamiltonian)
	endtime = time.time()
	return endtime - starttime
