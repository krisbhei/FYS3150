import numpy as np
import matplotlib.pyplot as plt

def plotLikelyState(filename):
    with open(filename,'r') as infile:
        line = infile.readline().split()
        temp = line[0]; trials = float(line[1])

        length = int(trials/100)
        energy_non_random = np.zeros(length) #because we sample for each 100 cycle
        energy_random = np.zeros(length)

        magnetization_non_random = np.zeros(length)
        magnetization_random = np.zeros(length)

        accepted_non_random = np.zeros(length)
        accepted_random = np.zeros(length)
        i = 0
        for line in infile:
            values_non_random = line.split()
            energy_non_random[i] = values_non_random[0]
            magnetization_non_random[i] = values_non_random[1]
            accepted_non_random[i] = values_non_random[2]

            next(infile)
            values_random = line.split()
            energy_non_random[i] = values_random[0]
            magnetization_non_random[i] = values_random[1]
            accepted_non_random[i] = values_random[2]
            i+=1
    plt.plot(accepted_non_random,energy_non_random,'ro')

if __name__ == "__main__":
    plotLikelyState("mostLikelyState_trials=100000_temp=1.dat")
