import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import os

def getAllResults():
    files = []
    for file in os.listdir("pos_vel"):
        if(file.endswith(".csv")):
            files.append(file)
    return files

def getMF_DM(file):

    split = file.split("_")
    print(split)
    return int(split[1]), float(split[3])
#
# for file in getAllResults():
#
#     mf, dm = getMF_DM(file)
#
#     data = pd.read_csv("pos_vel/" + file)
#
#     pos = data['pos'].values[::100]
#     vel = data['vel'].values[::100]
#     x = np.arange(0, len(pos) * 100)
#     plt.plot(pos, label="pos")
#     plt.plot(vel, label="vel")
#
#     if(mf > 0):
#         for vline in range(mf, len(x), mf):
#             plt.vlines(vline, -10, 10)
#     plt.show()
#


file = getAllResults()[4]

mf, dm = getMF_DM(file)

data = pd.read_csv("pos_vel/" + file)

pos = data['pos'].values[::100]
vel = data['vel'].values[::100]
x = np.arange(0, len(pos) * 100, 100)
plt.plot(x, pos, label="pos")
plt.plot(x, vel, label="vel")
plt.legend()
if (mf > 0):
     for vline in range(mf, len(x)*100, mf):
        plt.plot((vline,vline), (-10, 10), color='red')
plt.show()

#data = pd.read_csv("pos_vel_test.csv")

#
# #print(data.columns)
# #
# plt.figure(figsize=(26,9))
#
#
#
# pos = data['pos'].values[::100]
# vel = data['vel'].values[::100]
#
# x = np.arange(0, len(pos))
# plt.plot(x,pos, label="pos")
# plt.plot(x,vel, label="vel")
#
# plt.legend()
# plt.show()