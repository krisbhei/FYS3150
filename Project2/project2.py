#Wrapper file for calling the functions and generating the
#file discussed in the report.

import ctypes as ct
import numpy as np
from numpy.ctypeslib import ndpointer

from readPlotFile import read_file
from eig import find_eig
def run_non_time(rho_m1):
    n_values = [75,175,300,400,650]
    with open('non_compare_time.dat','w') as infile:
        for n in n_values:
            print "Taking time for solver with N = %d:"%n

            lib.non_interacting_time.argtypes=(ct.c_double,ct.c_int)
            lib.non_interacting_time.restype = (ct.c_double);

            time_c = lib.non_interacting_time(rho_max1,ct.c_int(n))
            time_p = find_eig(0,rho_m1,0,n)

            infile.write("\nTime used for non-interacting with N = %d:\n\
            - solver from jacobi.cpp used %.4f seconds\n\
            - numpy.linalg.eig from eig.py used %.4f seconds\n"\
            %(n,time_c,time_p))

            print "Done."

def run_non_interacting(rho_m1):
    rho_max1 = ct.c_double(rho_m1)

    n_values = [25,75,175,300,400]

    lib.non_interacting.argtypes=(ct.c_double,ct.c_int)
    lib.non_interacting.restype = ndpointer(dtype=ct.c_double,shape=(3,));
    filename = 'compare_non_rho=%.2f.dat'%rho_m1
    with open(filename,'w') as infile:
        infile.write("rho_max = %.2f\n"%rho_m1)
        for n in n_values:
            n_ = ct.c_int(n)

            res = lib.non_interacting(rho_max1,n_)
            infile.write("%d & %.6f & %.6f & %.6f \\\ \hline \n"%(n,res[0],res[1],res[2]))


def run_interacting(rho_m1,rho_m,coulomb=True):
    rho_max1 = ct.c_double(rho_m1)
    rho_max = ct.c_double(rho_m)
    filename = 'jacobi_lowest'
    if not coulomb:
        func_call = lib.interacting_case_noCoul
        filename +='_noCoul'
    else:
        func_call = lib.interacting_case
    filename += '.dat'
    func_call.argtypes = (ct.c_double,ct.c_double,ct.POINTER(ct.c_double),ct.c_int)

    w = [0.01,.5,1.,5.]
    c_arr = ct.c_double*len(w)
    w = c_arr(*w)

    rho_max = ct.c_double(50.)

    n = ct.c_int(300)
    w_n = ct.c_int(len(w))
    func_call(rho_max,rho_max1,w,n,w_n)

    print "Plotting the results"
    read_file(filename,coulomb)

def run_unit_tests():
    doublepp = np.ctypeslib.ndpointer(dtype=np.uintp)

    #test: finding largest absolute value of the e1onal elements
    #      in a symmetric matrix.
    dim = 7
    test_mat = np.zeros(shape=(dim,dim))
    max_expected = 6.6260
    e1 = [3.1415,-6.6260,1.6180,-2.7182,1.337,-0.42]
    for i in range(dim):
        for j in range(i):
            test_mat[i][j] = e1[j-1]
            test_mat[j][i] = test_mat[i][j]
    test_mat[3][3] = 15.
    test_mat[dim-1][dim-1] = 10.
    test_mat_pp = (test_mat.ctypes.data + np.arange(test_mat.shape[0]) * test_mat.strides[0])\
                 .astype(np.uintp)


    lib.test_find_max.argtypes = (doublepp,ct.c_int)
    lib.test_find_max.restype = (ct.c_double)
    max_comp = lib.test_find_max(test_mat_pp,ct.c_int(dim))
    assert abs(max_expected - max_comp)<1E-14, "Error in test_find_max: cant't find proper max value\nExpected: %.4f, got: %.4f"\
            %(max_expected,max_comp)

    #test: find the lowest 'eigenvalues' of a marix.
    #      will test on a random diagonal matrix
    dim = 10
    test_ind = np.eye(dim)
    diag = [5,6,2,9,1,8,3,8,10,-1]
    for i in range(dim):
        test_ind[i][i] =diag[i]

    test_ind_pp = (test_ind.ctypes.data + np.arange(test_ind.shape[0]) * test_ind.strides[0])\
                 .astype(np.uintp)

    test_expected_p = np.array([9,4,2]).astype(np.int32)

    lib.test_find_lowest.argtypes = (doublepp,ndpointer(ct.c_int),ct.c_int)
    lib.test_find_lowest.restype = (ct.c_int)

    res = lib.test_find_lowest(test_ind_pp,test_expected_p,ct.c_int(dim))
    assert res == 0, "Error in test_find_lowest: can't find indicies for the three lowest eigenvalues"

    #test: find eigenvalues and eigenvectors for a 3x3 matrix.
    #      also see if the eigenvectors are normalized and are orthogonal
    dim = 3
    test_eigvec = np.eye(dim)

    test_mat = np.zeros(shape=(dim,dim))
    test_mat[0] = [5,-2,0]
    test_mat[1] = [-2,6,2]
    test_mat[2] = [0,2,7]

    eig1 = 3; eig2 = 6;eig3 = 9
    eigvals = np.array([eig1,eig2,eig3]).astype(np.float64)

    e1 = 2./3
    e2 = -1./3
    test_eigvec[0] = [e1,e1,e2]
    test_eigvec[1] = [e1,e2,e1]
    test_eigvec[2] = [e2,e1,e1]

    test_mat_pp = (test_mat.ctypes.data + np.arange(test_mat.shape[0]) * test_mat.strides[0])\
                 .astype(np.uintp)

    test_eigvec_pp = (test_eigvec.ctypes.data + np.arange(test_eigvec.shape[0]) * test_eigvec.strides[0])\
                 .astype(np.uintp)

    lib.test_solver.argtypes = (doublepp,doublepp,ndpointer(ct.c_double),ct.c_int)
    lib.test_solver.restype = (ct.c_int)

    res = lib.test_solver(test_mat_pp,test_eigvec_pp,eigvals,dim)
    assert res != -3, "\nThe computed eigenvalues doesn't equal the expected eigenvalues";
    assert res != -2, "\nThe computed eigenvectors doesn't equal the expected eigenvectors";
    assert res != -1, "\nThe computed eigenvectors aren't normalized";

if __name__ == '__main__':
    lib = ct.cdll.LoadLibrary('./libjacobi.so')

    #Run unittests
    #--------------------------
    run_unit_tests()

    rho_max1 = 5.
    rho_max = 50.
    #Generate different values of the lowest eigvalues for different n
    #--------------------------

    #rhos = [1.,3.5,5,6.5]
    #for rho in rhos:
    #  run_non_interacting(rho)

    #run_non_time(rho_max1)

    #Generate files to plot for the interacting cases
    #--------------------------
    #run_interacting(rho_max1,rho_max)
    #run_interacting(rho_max1,rho_max,False)
