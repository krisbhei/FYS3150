import numpy as np
import matplotlib.pyplot as plt
import shutil,os

def plotLikelyState(x,y,temp,trials,quantifier,conf):

    plt.ylabel('$\\langle \\mathscr{%s} \\rangle$'%quantifier,fontsize=14)
    plt.xlabel('Number of Monte Carlo cycles')
    plt.ylim([np.min(y)-.001,np.max(y)+.001])
    plt.title("Finding the most likely state for an initial %s 20$\\times$20spin configuration with \n temperature %s and %d Monte Carlo cycles"%(conf,temp,trials))
    plt.plot(x,y,'k')
    plt.gca().grid(True)

    #plt.show()
    plt.savefig('plot_LikelyState_%s_%s_temp=%s_MC=%d.pdf'%(quantifier,conf,temp,trials))
    plt.clf()


def readLikelyState(filename):
    with open(filename,'r') as infile:
        line = infile.readline().split()
        temp = line[0]; trials = float(line[1])

        length = int(trials/100) #because we sample for each 100 cycle
        energy_non_random = np.zeros(length)
        energy_random = np.zeros(length)

        magnetization_non_random = np.zeros(length)
        magnetization_random = np.zeros(length)

        accepted_non_random = np.zeros(length)
        accepted_random = np.zeros(length)
        cycles = np.zeros(length)
        i = 0
        for line in infile:
            values_non_random = line.split()
            energy_non_random[i] = values_non_random[0]
            magnetization_non_random[i] = values_non_random[1]
            cycles[i] = values_non_random[2]
            accepted_non_random[i] = values_non_random[3]


            values_random = next(infile).split()
            energy_random[i] = values_random[0]
            magnetization_random[i] = values_random[1]
            accepted_random[i] = values_random[3]
            i+=1
    #Start: Plot the energies and magnetization
    plotLikelyState(cycles,energy_non_random,temp,trials,'E','ordered')
    plotLikelyState(cycles,magnetization_non_random,temp,trials,'|M|','ordered')
    plotLikelyState(cycles,energy_random,temp,trials,'E','random')
    plotLikelyState(cycles,magnetization_random,temp,trials,'|M|','random')
    #End: Plot the energies and magnetization

    #Start: Plot the accepted number of configurations per 100th step
    plt.plot(cycles,accepted_non_random,'ko',markersize=1)
    plt.title('Number of accepted configurations for \nan initial ordered spin configuration at temperature = %s '%temp)
    plt.ylabel('Percentage of accepted configurations per 100th cycle')
    plt.xlabel('Number of Monte Carlo cycles')
    plt.gca().grid(True)
    plt.savefig('plot_LikelyState_accepted_ordered_temp=%s_MC=%d.pdf'%(temp,trials))
    plt.clf()
    #plt.show()
    plt.plot(cycles,accepted_random,'ko',markersize=1)
    plt.title('Number of accepted configurations for \nan initial random spin configuration at temperature = %s '%temp)
    plt.ylabel('Percentage of accepted configurations per 100th cycle')
    plt.xlabel('Number of Monte Carlo cycles')
    plt.gca().grid(True)
    plt.savefig('plot_LikelyState_accepted_ordered_temp=%s_MC=%d.pdf'%(temp,trials))
    plt.clf()
    #plt.show()
    #End: Plot the accepted number of configurations per 100th step

if __name__ == "__main__":

    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    params = {'text.latex.preamble' : [r'\usepackage{mathrsfs}']}
    plt.rcParams.update(params)

    readLikelyState("mostLikelyState_trials=1000000_temp=1.dat")
    readLikelyState("mostLikelyState_trials=1000000_temp=2.4.dat")
