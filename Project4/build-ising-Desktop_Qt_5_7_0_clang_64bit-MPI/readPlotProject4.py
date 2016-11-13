import numpy as np
import matplotlib.pyplot as plt
import shutil,os
import numpy as np
import matplotlib.pyplot as plt
import shutil,os

def readSearchEnergy(filename,ordering):
    with open(filename,'r') as infile:
        temp = infile.readline().strip('\n')
        energies = [float(val) for val in infile.readline().split()]
        numEnergies = len(energies)
        hist = [float(infile.readline()) for i in xrange(numEnergies)]

        #We are at variance;
        varianceE = float(infile.readline())



    print varianceE
    plt.title('Number of occurence of possible energies for 20$\\times$20 lattice at \ntemperature T=%s and %s initial ordering'%(temp,ordering))
    plt.xlabel('Possible energies')
    plt.ylabel('Number of occurence')
    plt.plot(energies,hist,'k')
    plt.xlim([-800,-200])
    plt.savefig('plot_'+filename[:-4]+'.pdf')
    plt.show()

def plotPhase(temp,q,title,xlabel,ylabel,fig_num):
    plt.figure(fig_num)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel,fontsize=14)
    plt.plot(temp,q)
    plt.hold('on')


def readPhaseTransitions(filename):
    with open(filename,'r') as infile:
        temps = [float(value) for value in infile.readline().split()]
        num_temps = len(temps)
        labels = []
        for L in [40,60,100,140]:
            avgE = np.zeros(num_temps)
            avgAbsM = np.zeros(num_temps)
            Cv = np.zeros(num_temps)
            X = np.zeros(num_temps)

            for i in range(num_temps):
                values = [float(v) for v in infile.readline().split()]
                avgE[i] = values[0]
                Cv[i] = values[1]
                avgAbsM[i] = values[2]
                X[i] = values[3]
            print X
            labels.append('%d$\\times$%d'%(L,L))
            plotPhase(temps,avgE,'Mean energy for temperature in [2,2.3] for different L$\times$L lattices','Temperature','$\\langle \\mathscr{E} \\rangle$',1)
            plotPhase(temps,Cv,'Heat capacity for temperature in [2,2.3] for different L$\times$L lattices','Temperature','$\\mathscr{C}_V $',2)
            plotPhase(temps,avgAbsM,'Mean magnetization for temperature in [2,2.3] for different L$\times$L lattices','Temperature','$\\langle \\mathscr{|M|} \\rangle$',3)
            plotPhase(temps,X,'Susceptibility for temperature in [2,2.3] for different L$\times$L lattices','Temperature','$\\mathscr{X}$',4)
    print labels
    for i in xrange(4):
        plt.figure(i+1)
        ax = plt.subplot(111)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=5)
        plt.legend(labels)
        plt.savefig('plot_readPhase%d.pdf'%i)

def plotLikelyState(x,y,temp,trials,quantifier,conf):
    # plt.subplots_adjust(hspace=.4)
    # plt.subplot(2,1,1)
    # plt.ylabel('$\\langle \\mathscr{%s} \\rangle$'%quantifier,fontsize=14)
    # plt.xlabel('Number of Monte Carlo cycles')
    # plt.ylim([np.min(y)-.001,np.max(y)+.001])
    plt.title("Finding the most likely state for an initial %s 20$\\times$20 spin configuration with \n temperature %s and %d Monte Carlo cycles"%(conf,temp,trials))
    #
    # plt.plot(x,y)
    # plt.gca().grid(True)
    #
    # plt.xlim([0,650000])
    # plt.ylim([-2,-1.95])
    # #plt.ylim([.97,1])
    #
    # plt.subplot(2,1,2)
    plt.ylabel('$\\langle \\mathscr{%s} \\rangle$'%quantifier,fontsize=14)
    plt.xlabel('Number of Monte Carlo cycles')
    #plt.title("Zoom of the plot above")

    #plt.ylim([.9990,.9995])
    #plt.ylim([-1.9975,-1.9965])

    #plt.xlim([0,650000])
    plt.gca().grid(True)
    plt.plot(x,y)
    #plt.show()
    plt.savefig('zoomx_0_zoomy_0_plot_LikelyState_%s_%s_temp=%s_MC=%d.pdf'%(quantifier,conf,temp,trials))
    plt.clf()


def readLikelyState(filename,ylim=None):
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
    #plotLikelyState(cycles,energy_random,temp,trials,'E','random')
    #plotLikelyState(cycles,magnetization_random,temp,trials,'|M|','random')
    #End: Plot the energies and magnetization

    #Start: Plot the accepted number of configurations per 100th step
    # plt.plot(cycles,accepted_non_random,'k')
    # plt.title('Number of accepted configurations for \nan initial ordered spin configuration at temperature = %s '%temp)
    # plt.ylabel('Percentage of accepted configurations per 100th cycle')
    # plt.xlabel('Number of Monte Carlo cycles')
    # plt.xlim([0,800000])
    # plt.gca().grid(True)
    # plt.savefig('zoom_800000_plot_LikelyState_accepted_ordered_temp=%s_MC=%d.pdf'%(temp,trials))
    # plt.show()
    # plt.clf()

    # plt.plot(cycles,accepted_random)
    # plt.title('Number of accepted configurations for \nan initial random spin configuration at temperature = %s '%temp)
    # plt.ylabel('Percentage of accepted configurations per 100th cycle')
    # plt.xlabel('Number of Monte Carlo cycles')
    # plt.gca().grid(True)
    # plt.xlim([0,800000])
    # plt.ylim([0.07,0.105])
    # plt.savefig('xzoom_800000_yzoom007+0105_plot_LikelyState_accepted_random_temp=%s_MC=%d.pdf'%(temp,trials))
    # plt.show()
    # plt.clf()
    #End: Plot the accepted number of configurations per 100th step
def plotAndReadForLikelyState():
    #readLikelyState("mostLikelyState_trials=1000000_temp=1.dat")
    readLikelyState("mostLikelyState_trials=10000000_temp=2.4.dat")

def plotAndReadSearchEnergy():
    temps = [1,2.4]
    orderings = ['ordered','random']
    base = 'searchEnergies_temp='
    for temp in temps:
        for init in orderings:
            filename = base + '%.6f'%(temp) + '_%s'%(init)+'.dat'
            readSearchEnergy(filename,init)
if __name__ == "__main__":

    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    params = {'text.latex.preamble' : [r'\usepackage{mathrsfs}']}
    plt.rcParams.update(params)
    plotAndReadForLikelyState()
    #plotAndReadSearchEnergy()
    #readPhaseTransitions('phaseTransitions_Tstart=2_Tend=2.3_Tstep=0.05.dat')
