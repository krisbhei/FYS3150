#Wrapper file which calls the functions discussed in the report.

import ctypes as ct
import numpy as np
from numpy.ctypeslib import ndpointer
from matplotlib import pyplot as plt

from solarsim import *

C_DOUBLE = ct.c_double; C_CHAR = ct.c_char_p; C_INT = ct.c_int

def run_relativist_correction(start_yr,end_yr,filename,n,savefigure=False):
    s_yr = C_DOUBLE(start_yr); e_yr = C_DOUBLE(end_yr)
    c_filename = filename.encode('utf-8')

    c_n = C_INT(n)
    generalRelativity = lib.generalRelativity
    generalRelativity.argtypes=(C_DOUBLE,C_DOUBLE,C_CHAR,C_INT)
    generalRelativity.restype=(None)

    generalRelativity(s_yr,e_yr,c_filename,c_n)
    time = []
    perihelions = []
    totalEnergy = []
    with open(filename,'r') as infile:
        line = infile.readline().split()
        while(line != []):
            time.append(float(line[-1]))
            perihelions.append(float(infile.readline().split()[-1]))
            totalEnergy.append(float(infile.readline().split()[-1]))
            line = infile.readline().split()

    plt.figure(1)
    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    plt.title('Perihelions in arcseconds found plotted against \nthe times they were found during %d year(s)'%end_yr)
    plt.plot(time,perihelions)
    plt.xlabel('Time')
    plt.ylabel(r'$\theta_p$',fontsize=18)

    plt.figure(2)
    plt.title('Total energy for every found perihelion')
    plt.plot(time,totalEnergy)
    plt.xlabel('Time')
    plt.ylabel(r'$E_{kinetic} + E_{potential}$',fontsize=14)

    if savefigure:
        plt.figure(1)
        plt.savefig('n=%d_perihelion_relativistic.pdf'%n)
        plt.figure(2)
        plt.savefig('n=%d_totenergy_relativistic.pdf'%n)
    else:
        plt.show()

    return #End: run_relativist_correction
def run_sun_earth_jupiter(inc_factor,start_yr,end_yr,filename,n,savefigure=False):
    s_yr = C_DOUBLE(start_yr); e_yr = C_DOUBLE(end_yr)
    c_inc_f = C_DOUBLE(inc_factor)
    c_n = C_INT(n)

    c_filename = filename.encode('utf-8')

    sunEarthJupiter = lib.sunEarthJupiter
    sunEarthJupiter.argtypes=(C_DOUBLE,C_DOUBLE,C_DOUBLE,C_CHAR,C_INT)
    sunEarthJupiter.restype=(None)
    sunEarthJupiter(c_inc_f,s_yr,e_yr,c_filename,c_n)

    plotPosition(filename,'Three body system with Earth, Jupiter and Sun using N = %d.'%n + \
                '\nIncreasing factor of Jupiter\'s mass: '+str(inc_factor), \
                savefigure)

    t = np.linspace(start_yr,end_yr,n)
    k,p,a,m = getConservationEnergyAndMomentum(filename)

    tot = k+p
    at = [euclideanNorm(m) for m in a]
    plotEnergyAndMomentum(t,at,tot,'Verlet for %d years and increasing factor: %d'%(end_yr,inc_factor),savefigure)

    return #End: run_sun_earth_jupiter

def run_escape_sun(vx,vy,start_yr,end_yr,filename,n,savefigure=False):
    s_yr = C_DOUBLE(start_yr); e_yr = C_DOUBLE(end_yr)
    c_n = C_INT(n)

    c_filename = filename.encode('utf-8')

    escapeSun = lib.escapeSun
    escapeSun.argtypes=(C_DOUBLE,C_DOUBLE,C_DOUBLE,C_DOUBLE,C_CHAR,C_INT)
    escapeSun.restype=(None)

    c_vx = C_DOUBLE(vx)
    c_vy = C_DOUBLE(vy)

    escapeSun(c_vx,c_vy,s_yr,e_yr,c_filename,c_n)

    filename_experiment = 'velocity_experiment_'+filename;
    plotPosition(filename_experiment,'Finding inital velocity for Planet to escape Sun:\n $v_{x}$ = '+str(vx)+\
                ' $v_{y}$ = '+str(vy),savefigure)

    filename_exact = 'velocity_exact_'+filename
    plotPosition(filename_exact,'Using exact inital velocity for Planet to escape Sun:\n $v_{x}$ = 0,'+\
                ' $v_{y}$ = $\sqrt{2GM_{\odot}/r} = \sqrt{8\pi^2}$',savefigure)
    return #End: run_escape_sun

def run_earth_sun(start_yr,end_yr,filename,n,savefigure=False):
    s_yr = C_DOUBLE(start_yr); e_yr = C_DOUBLE(end_yr)
    c_n = C_INT(n)

    c_filename = filename.encode('utf-8')

    testEarthSun = lib.testEarthSun
    testEarthSun.argtypes=(C_DOUBLE,C_DOUBLE,C_CHAR,C_INT)
    testEarthSun.restype=(None)
    testEarthSun(s_yr,e_yr,c_filename,c_n)

    #Now we have two files; one for Euler and one for Verlet
    t = np.linspace(start_yr,end_yr,n)

    file_e = 'e_'+filename
    file_v = 'v_'+filename
    files = [file_e,file_v]
    odes = ['Euler','Verlet']
    for i in xrange(2):
        print "-Using method: ",odes[i]
        plotPosition(files[i],('Solved orbit of Earth using %s \nwith N = '+str(n))%odes[i], \
                    savefigure)

        k,p,a,m = getConservationEnergyAndMomentum(files[i])
        tot = k+p
        at = [euclideanNorm(m) for m in a]
        plotEnergyAndMomentum(t,at,tot,odes[i],savefigure)

        testConservationEnergyAndMomentum(k,p,a)


def run_all_planets(start_yr,end_yr_all,end_yr_three,filename,n_all,n_three,savefigure=False):


    s_yr = C_DOUBLE(start_yr); e_yr_all = C_DOUBLE(end_yr_all);e_yr_three = C_DOUBLE(end_yr_three)
    c_filename = filename.encode('utf-8')

    c_n_all = C_INT(n_all)
    c_n_three = C_INT(n_three)


    allPlanets = lib.allPlanets
    allPlanets.argtypes=(C_DOUBLE,C_DOUBLE,C_DOUBLE,C_CHAR,C_INT,C_INT)
    allPlanets.restype=(None)
    allPlanets(s_yr,e_yr_all,e_yr_three,c_filename,c_n_all,c_n_three)

    #Start: Plot Jupiter - Sun - Earth and test for zero momentum
    threeBody_file = 'threeBodyZeroMomentum_'+filename

    #plotPosition(threeBody_file,'Earth, Jupiter and Sun during %d year(s) and N = %d.'%(end_yr_three,n_three) + \
    #              '\nThe Sun has velocity such that the total momentum of the system is zero', savefigure)

    momentum = getConservationEnergyAndMomentum(threeBody_file)[-1]
    testForZeroMomentum(momentum)
    #End: Plot Jupiter - Sun - Earth and test for zero momentum

    #Plot system for all planets
    plotPosition(filename, 'All planets in the solar system after %d year(s)\n with N = %d'%(end_yr_all,n_all),savefigure)

    return #End: run_all_planets

def run_time_euler_verlet(start_yr,end_yr,filename,n,num_time_taking):
    s_yr = C_DOUBLE(start_yr); e_yr = C_DOUBLE(end_yr)
    c_n = C_INT(n)

    timeEuler = lib.timeEuler
    timeEuler.argtypes=(C_DOUBLE,C_DOUBLE,C_INT)
    timeEuler.restype=(C_DOUBLE)

    timeVerlet = lib.timeVerlet
    timeVerlet.argtypes=(C_DOUBLE,C_DOUBLE,C_INT)
    timeVerlet.restype=(C_DOUBLE)

    time_e = np.zeros(num_time_taking)
    time_v = np.zeros(num_time_taking)
    average_e = 0
    average_v = 0
    for i in xrange(num_time_taking):
            time_e[i] = timeEuler(s_yr,e_yr,c_n)
            time_v[i] = timeVerlet(s_yr,e_yr,c_n)
            average_e += time_e[i]
            average_v += time_v[i]

    avarage_e = average_e/num_time_taking
    avarage_v = average_v/num_time_taking
    with open('n=%d_time_euler_vs_verlet.dat'%n,'w') as outfile:
        outfile.write('Time comparison in seconds between Euler and Verlet for N = %d\n'%n)
        outfile.write('----------------------------------------\n')
        for i in xrange(num_time_taking):
            outfile.write('Time for %d-th call: \n'%i)
            outfile.write('           - Euler: %.8f\n'%time_e[i])
            outfile.write('           - Verlet: %.8f\n'%time_v[i])
        outfile.write('----------------------------------------\n')
        outfile.write('Average time: \n')
        outfile.write('Euler: %.8f seconds\n'%avarage_e)
        outfile.write('Verlet: %.8f seconds'%avarage_v)

    return #End: run_time_euler_verlet

if __name__ == '__main__':
    lib = ct.cdll.LoadLibrary('./libsolarsim.so')
    run = True
    while run:
        print "[1] testEarthSun:          run Euler and Verlet for 10^2, 10^4 and 10^5 meshpoints and check for conservation of energy and momentum\n"+\
        "[2] timeEulerVerlet:       runs the Earth-Sun system 5 times for Euler and Verlet and takes time for each call. \n"+ \
        "[3] escapeSun:             runs the Earth-Sun system by giving Earth an experimental escape velocity and exact escape velocity\n"+\
        "[4] sunEarthJupiter:       runs the Earth-Sun-Jupiter system for 10^2,10^3 and 10^4 meshpoints and increases the mass of Jupiter by 1, 10^1 and 10^3 for every meshpoint \n" +\
        "[5] allPlanets:            runs Sun-Earth_Jupiter by giving the Sun an initial velocity, then runs for all planets with inital data from NASA\n"+ \
        "[6] generalRelativity:     runs the Sun-Mercury system with relativistic correction"+\
        "\n[anything else] exit"

        try:
            arg=input("Call the desired function by typing its corresponding number: ")

            if arg == 1:
                #Testing the Earth-Sun system:
                n_exp = [2,4,5]
                for n in n_exp:
                    points = 10**n
                    print "\nRunning Earth-Sun system for N = %d for 1 year"%points
                    filename = ('n=%d_')%points + 'earth-sun.dat'
                    run_earth_sun(0,1,filename,points,savefigure=False)

                print "- Done running testEarthSun -"

            elif arg == 2:
                #Taking time of Euler and Verlet:
                n_exp = [6,7,8]
                number_iter = 5
                for n in n_exp:

                    points = 10**n
                    print "Taking time for Euler and Verlet using N = 10^%d meshpoints" %n

                    filename = 'n=%d_time.dat'%points
                    run_time_euler_verlet(0,1,filename,points,number_iter,savefigure=False)

                print "- Done timeEulerVerlet -"

            elif arg == 3:
                #Testing the speed for a planet to escape the Sun:
                #(0,2.9pi,0) is the velocity found by experimenting by which causes a
                #planet with 1 AU from the Sun to escape.
                run_escape_sun(0,2.9*np.pi,0,1,"escapeSun.dat",10**4,savefigure=False)

                print "- Done escapeSun -"


            elif arg == 4:
                #Test for how the increasing mass of Jupiter will affect
                #the orbit of the Earth
                filename = 'sun_earth_jupiter_incF='
                n_exp = [2,3,4]
                years = 15
                print 'Running for %d years'%years
                for n in n_exp:
                    N = 10**n
                    for i in [0,1,3]:
                        inc = 10**i
                        print 'Running the Sun-Earth-Jupiter system using N = %d and increasing Jupiter\'s mass by %d'%(N,inc)
                        f = ('n=%d'%N)+filename+str(inc)
                        run_sun_earth_jupiter(inc,0,years,f+'.dat',N,savefigure=False)

                N = 10**4
                print "- Done sunEarthJupiter -"

            elif arg == 5:
                #First, run the simulation for three bodies; Earth, Jupiter and Sun with
                #the Sun having an initial velocity.
                #Then, run the simulation for all 8 planets in the solar system
                #including Pluto with initial velocity and positon from NASA
                n = 6
                n1 = 1
                years_all = 100
                years_three = 1

                N_all = 10**n
                N_three = 10**n1
                print 'Running three-body system with N = %d for %d year(s)'%(N_three,years_three)+ \
                      '\nRunning all the planets with N = %d for %d year(s)'%(N_all,years_all)
                filename = 'yr=%d_n=%d_'%(years_all,N_all) + 'allPlanets.dat'

                run_all_planets(0,years_all,years_three,filename,N_all,N_three,savefigure=False)
                print '- Done allPlanets -'

            elif arg == 6:
                #Run the simulation with Sun and Mercury to see how the relativistic correction
                #makes the orbit of Mercury elliptic in contrast to using pure Newtonian gravitational force
                run_relativist_correction(0,100,'RelativisticCorrection.dat',10**8,savefigure=False)
            else:
                run = False
        except ValueError:
            run = False
        except SyntaxError:
            run = False
