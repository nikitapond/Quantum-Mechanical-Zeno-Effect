import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("results/average_amplitudes2.csv")



def plotAllFreq():
    '''

    Plots each amplitude against max disturbance for different measurement frequencies

    '''

    freq = data['MFreqs'].values
    plt.figure(figsize=(12,9))
    x = data.columns.values[1:]

    for mf in freq:
        conFreq = data.loc[data['MFreqs'] == mf].values[0][1:]
        plt.plot(x, conFreq, label="inf" if mf==-1 else mf)
    plt.xlabel("Max Disturbance (m/s)")
    plt.ylabel("Average Amplitude (m)")
    plt.legend(title="Time steps between disturbances")
    plt.savefig("results/amp_against_disturbance.pdf")
    plt.show()
plotAllFreq()

def plotAllDist():

    mDist = data.columns.values[1:]
    plt.figure(figsize=(12,9))

    xticks = data['MFreqs'].values
    xticks = np.where(xticks == -1, "Inf", xticks)

    for d in mDist[::2]:
        vals = data[d].values
        plt.plot(vals, label=d)
    plt.xticks(np.arange(0, 8), labels= xticks)

    plt.legend(title="Max Disturbance (m/s)")
    plt.ylabel("Average Amplitude (m)")
    plt.xlabel("Time steps between disturbances")
    plt.savefig("results/amp_against_dist_freq.pdf")
    plt.show()

    print(mDist)

plotAllDist()



#
# def plotConstFreq(freq):
#
#     conFreq = data.loc[data['MFreqs'] == freq].values[0][1:]
#     print(conFreq)
#     plt.plot(conFreq)
#
#
#
# plotConstFreq(200)
# plt.show()