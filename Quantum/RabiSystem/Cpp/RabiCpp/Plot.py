import numpy as np
import matplotlib.pyplot as plt

FONT_SIZE = 14
fig = plt.figure(figsize=(12,9))
ax = fig.subplots()



measure_f = [-1, 0.2, 0.1, 0.04,16/500, 8/500, 4/500, 2/500]
ax.set_title("Plot of average MCFW results over average of 20000 iterations \nfor different measurement frequencies", fontsize=FONT_SIZE)

time = np.linspace(0, 4, 2000)
for mf in measure_f:

    label = "inf" if mf<0 else "%.3f" % round(mf, 3)
    file = "mcwf_rabi_res_mf_" + label + ".csv"
    print(file)
    data = np.genfromtxt(file, delimiter=',')
    plt.plot(time, data[1], label = "$\Delta M_f=$" + label)

ax.set_ylabel("Probability of 'b' state occupancy", fontsize=FONT_SIZE)
ax.set_xlabel(r"Time (units of Rabi Frequency $R_{Hz}$)", fontsize=FONT_SIZE)
plt.legend(loc='upper left',  fontsize=FONT_SIZE)


plt.savefig("mcfw_high_it_results.png")
plt.savefig("mcfw_high_it_results.pdf")


plt.show()
# print(data.shape)
# plt.plot(data[0])
# plt.plot(data[1])
#
# plt.show()