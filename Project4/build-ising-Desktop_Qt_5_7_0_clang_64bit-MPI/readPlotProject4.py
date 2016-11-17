import numpy as np
import matplotlib.pyplot as plt
import shutil,os
import numpy as np
import matplotlib.pyplot as plt
import shutil,os

def findA(critTemp,Ls):
    #start: find the a's
    a = np.zeros(shape=(4,3))
    i = 0
    for temp in critTemp:
        for _i in xrange(3):
            j = (i+1+_i)%4
            dividend = 1./( (1./Ls[i])-(1./Ls[j]) )
            a[i,_i] = (temp - critTemp[j])*dividend
        i += 1
    #End: find the a's
    _a = np.sum(np.sum(a))/12.
    computedCritTemp = np.sum([critTemp[i] - _a/float(Ls[i]) for i in xrange(4)])/4.

    return _a,a,computedCritTemp

def critTemp():
    Ls = [40,60,100,140]
    critTemp = np.zeros(len(Ls))
    critI = 0
    plt.xlabel('Temperature')
    plt.ylabel('$\mathscr{C}_V$',fontsize=14)
    plt.xlim([2.18,2.42])
    plt.title('Heat capacity for temperatures in [2.1,2.5] for different sizes of lattices')

    labels = []
    critTempsSusIndex = 0
    maxSus = 0
    critTempsSus = np.zeros(len(Ls))
    for L in Ls:
        filename = 'new_criticalTemp_Tstart=2.1_Tend=2.5_%d.dat'%L
        with open(filename,'r') as infile:
            temps = [float(T) for T in infile.readline().split()]
            T_len = len(temps)
            heatCapasity = np.zeros(T_len)
            susceptibility = np.zeros(T_len)
            magnetization = np.zeros(T_len)
            infile.readline()

            for i in xrange(T_len):
                values = infile.readline().split()
                heatCapasity[i] = float(values[1])
                susceptibility[i] = float(values[3])
                if susceptibility[i] > maxSus:
                    maxSus = susceptibility[i]
                    critTempsSusIndex = i
                magnetization[i] = float(values[2])

            critTempsSus[critI] = temps[critTempsSusIndex]
            critTemp[critI] = float(infile.readline())
            critI += 1
        plt.plot(temps,heatCapasity)
        labels.append('$%d\\times%d$'%(L,L))
        plt.hold('on')
    plt.legend(labels)
    plt.savefig('1all_plot_'+filename[:-4]+'.pdf')
    # plt.clf()


    _a,a,_critTemp = findA(critTemp,Ls)

    print _critTemp
    print findA(critTempsSus,Ls)[-1]


    # xis = np.zeros(4)
    # for i in xrange(4):
    #     xis[i] = 1./abs(computedCritTemp - critTemp[i])
    #
    # __a = np.zeros(shape=(4,3))
    #
    # i = 0
    # for x in xis:
    #     for _i in xrange(3):
    #         j = (i+1+_i)%4
    #         dividend = 1./( (1./Ls[i])-(1./Ls[j]) )
    #         __a[i,_i] = (x - xis[j])*dividend
    #     i += 1
    #
    # newA = np.sum(np.sum(__a))/12.
    #
    # computedCritTemp2 = np.sum([critTemp[i] - newA/float(Ls[i]) for i in xrange(4)])/4.
    # print computedCritTemp2


    #Write the values of the critical temperatures and different a's in a table in latex format
    with open('1tableCritTemp.txt','w') as outfile:
        outfile.write('%.5f '%critTemp[0])
        for i in xrange(3):
            outfile.write('& %.5f '%critTemp[i+1])
        outfile.write('\\\\ \\hline\n\n')

        for i in xrange(4):
            outfile.write('$%d - L_j$ '%Ls[i])
            for j in xrange(3):
                outfile.write('& %.5f '%(a[i,j]))
            outfile.write('\\\\ \\hline\n')
        outfile.write('%.5f'%(_a))
    #End: write table in Latex format

    #Start: compute
def readSearchEnergy(filename,ordering,temp):
    with open(filename,'r') as infile:
        temp = infile.readline().strip('\n')
        energies = [float(val) for val in infile.readline().split()]
        numEnergies = len(energies)
        hist = np.array([float(infile.readline()) for i in xrange(numEnergies)])

        #We are at the standard deviation;
        stdE = float(infile.readline())

    hist /= np.sum(hist)
    plt.subplots_adjust(hspace=.4)
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.xlim([-2.1,2])
    plt.title('Number of occurence of possible energies for 20$\\times$20 lattice at \ntemperature T=%s and %s initial ordering'%(temp,ordering))
    plt.xlabel('Possible energies')
    plt.ylabel('Probability of occurence')
    markerline, stemlines, baseline = plt.stem(energies, hist,'-')
    plt.setp(markerline, 'markerfacecolor', 'k','markersize',2)
    plt.setp(baseline, 'color', 'k', 'linewidth', 2)
    plt.setp(stemlines,'color','k')

    plt.subplot(2,1,2)
    #plt.xlim([-2.1,-1.9])  #For low temperature
    plt.xlim([-1.7,-.9])   #For high temperature
    plt.title('Zoom of the plot above')
    plt.xlabel('Possible energies')
    plt.ylabel('Probability of occurence')
    markerline, stemlines, baseline = plt.stem(energies, hist,'-')
    plt.plot(energies,hist,'k--')
    plt.setp(markerline, 'markerfacecolor', 'k')
    plt.setp(baseline, 'color', 'k', 'linewidth', 2)
    plt.setp(stemlines,'color','k')

    plt.savefig('zoom_-2-1.9_plot_'+filename[:-4]+'.pdf')
    plt.show()

def plotPhase(temp,q,title,xlabel,ylabel,l,fig_num):
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

            labels.append('%d$\\times$%d'%(L,L))
            legend = '%d$\\times$%d'%(L,L)
            plotPhase(temps,avgE,'Mean energy for temperature in [2,2.3] for different L$\\times$L lattices','Temperature','$\\langle \\mathscr{E} \\rangle$',legend,0)
            plotPhase(temps,Cv,'Heat capacity for temperature in [2,2.3] for different L$\\times$L lattices','Temperature','$\\mathscr{C}_V $',legend,1)
            plotPhase(temps,avgAbsM,'Mean magnetization for temperature in [2,2.3] for different L$\\times$L lattices','Temperature','$\\langle \\mathscr{|M|} \\rangle$',legend,2)
            plotPhase(temps,X,'Susceptibility for temperature in [2,2.3] for different L$\\times$L lattices','Temperature','$\\mathscr{X}$',legend,3)
    plt.figure(0)
    plt.legend(labels,loc='upper left')
    plt.savefig('1plot_readPhase0.pdf')
    plt.figure(1)
    plt.legend(labels,loc='upper left')
    plt.savefig('1plot_readPhase1.pdf')
    plt.figure(2)
    plt.legend(labels)
    plt.savefig('1plot_readPhase2.pdf')
    plt.figure(3)
    plt.legend(labels,loc='upper left')
    plt.savefig('1plot_readPhase3.pdf')
    #plt.savefig('plot_readPhase%d.pdf'%i)
    #plt.show()
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
    #plt.ylim([-1.25,-1.23])
    #plt.xlim([0,650000])
    #plt.ylim([.45,.5])
    plt.ylim([.425,.475])
    plt.gca().grid(True)
    plt.plot(x,y)
    #plt.show()
    plt.savefig('zoomx_0_zoomy_.425+.475_plot_LikelyState_%s_%s_temp=%s_MC=%d.pdf'%(quantifier,conf,temp,trials))
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
    #plotLikelyState(cycles,energy_non_random,temp,trials,'E','ordered')
    #plotLikelyState(cycles,magnetization_non_random,temp,trials,'|M|','ordered')
    #plotLikelyState(cycles,energy_random,temp,trials,'E','random')
    #plotLikelyState(cycles,magnetization_random,temp,trials,'|M|','random')
    #End: Plot the energies and magnetization

    #Start: Plot the accepted number of configurations per 100th step
    plt.plot(cycles,accepted_non_random)
    plt.title('Accepted configurations for \nan initial ordered spin configuration at temperature = %s '%temp)
    plt.ylabel('Percentage of accepted configurations per 100th cycle')
    plt.xlabel('Number of Monte Carlo cycles')
    plt.xlim([0,800000])
    #plt.ylim([26,27])
    plt.ylim([26.75,27.25])
    plt.gca().grid(True)
    plt.savefig('zoomy_26.75+27.25_plot_LikelyState_accepted_ordered_temp=%s_MC=%d.pdf'%(temp,trials))
    # plt.show()
    # plt.clf()
    #
    # plt.plot(cycles,accepted_random)
    # plt.title('Accepted configurations for \nan initial random spin configuration at temperature = %s '%temp)
    # plt.ylabel('Percentage of accepted configurations per 100th cycle')
    # plt.xlabel('Number of Monte Carlo cycles')
    # plt.gca().grid(True)
    # plt.xlim([0,800000])
    # #plt.ylim([0.07,0.105])
    # #plt.ylim([26.5,27.2])
    # plt.ylim([0,30])

    plt.savefig('zoomy_30_plot_LikelyState_accepted_random_temp=%s_MC=%d.pdf'%(temp,trials))
    # plt.show()
    # plt.clf()
    #End: Plot the accepted number of configurations per 100th step
def plotAndReadForLikelyState():
    readLikelyState("mostLikelyState_trials=1000000_temp=1.dat")
    readLikelyState("lmostLikelyState_trials=10000000_temp=2.4.dat")

def plotAndReadSearchEnergy():
    temps = [1,2.4]
    orderings = ['ordered']
    base = 'searchEnergies_temp='
    for temp in temps:
        for init in orderings:
            filename = base + '%.6f'%(temp) + '_%s'%(init)+'.dat'
            readSearchEnergy(filename,init,temp)


if __name__ == "__main__":

    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    params = {'text.latex.preamble' : [r'\usepackage{mathrsfs}']}
    plt.rcParams.update(params)
    #plotAndReadForLikelyState()
    #plotAndReadSearchEnergy()
    #readPhaseTransitions('phaseTransitions_Tstart=2_Tend=2.4_Tstep=0.015.dat')
    critTemp()
