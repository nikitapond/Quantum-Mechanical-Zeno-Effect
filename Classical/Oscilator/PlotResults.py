import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("all_results.csv")


def plotAllFreq():

    freq = data['MFreqs'].values
    plt.figure(figsize=(12,9))
    for mf in freq:
        conFreq = data.loc[data['MFreqs'] == mf].values[0][1:]
        plt.plot(conFreq, label=mf)
    plt.legend()
    plt.show()


plotAllFreq()
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