#Module to plot the values generated in a file.
#Coded to read files created by extern "C" void interacting_case in jacobi.cpp.
#Called in project2.py


from matplotlib import pyplot as plt
import numpy as np

def plot_values(rho,u,w_r,lamb,n,k,coulomb=True):
    plt.figure(k)
    plt.plot(rho,u)
    plt.xlabel(r'$\rho$',fontsize=16)
    plt.ylabel(r'probability density $u^2$',fontsize=14)
    title = 'Interacting case for N = %d using the \
            \n eigenstate with smallest energylevel $\lambda = $%s and $\omega_r = $%s'\
            %(n,lamb,w_r)

    filename = 'interacting_wr=%s_n=%d'%(w_r,n)

    if not coulomb:
        title += '\nwithout Coulomb interaction'
        filename += '_nCoul'
    filename+='.png'

    plt.title(title,fontsize=14)
    plt.savefig(filename)
    plt.clf()

def read_file(filename,coulomb=True):
    with open(filename,'r') as infile:
            plt.rc('text',usetex=True)
            plt.rc('font', family='serif')

            n = int(infile.readline().split()[-1])

            w_r = infile.readline().split()[-1]
            lamb = infile.readline().split()[-1]
            values = infile.readline().split()

            rho = np.zeros(n)
            u = np.zeros(n)
            i = 0
            k = 0 #number of plots
            while(values != []):
                try:
                    rho[i] = float(values[0])
                    u[i] = float(values[1])
                    i+=1
                    values = infile.readline().split()
                except ValueError: #we have reached new values for wr and possible rho
                    plot_values(rho,u,w_r,lamb,n,k,coulomb)
                    k+=1
                    i = 0
                    w_r = values[-1]
                    lamb = infile.readline().split()[-1]
                    values = infile.readline().split()

            plot_values(rho,u,w_r,lamb,n,k,coulomb)

if __name__ == '__main__':
    read_file('jacobi_lowest.dat')
