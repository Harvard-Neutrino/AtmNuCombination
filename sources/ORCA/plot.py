import matplotlib.pyplot as plt
import numpy as np

import util

def plot_zenith_errors():
    e, ebar, mu, mubar = util.get_zenith_error()

    x = np.linspace(1.8, 53, 100000)
    ye = np.zeros_like(x)
    yeb = np.zeros_like(x)
    ymu = np.zeros_like(x)
    ymub = np.zeros_like(x)

    for i in range(len(x)):
        ye[i] = e(x[i])
        yeb[i] = ebar(x[i])
        ymu[i] = mu(x[i])
        ymub[i] = mubar(x[i])

    plt.plot(x, ye, color = "red", linewidth = 1.8, linestyle = '-', label = "Electron Neutrino")
    plt.plot(x, yeb, color = "red", linewidth = 1.8, linestyle = '--', label = "Electron AntiNeutrino")
    plt.plot(x, ymu, color = "blue", linewidth = 1.8, linestyle = '-', label = "Muon Neutrino")
    plt.plot(x, ymub, color = "blue", linewidth = 1.8, linestyle = '--', label = "Muon AntiNeutrino")

    plt.xscale('log')

    plt.xlabel("Neutrino Energy [Gev]")
    plt.ylabel("Median Angular Error [deg]")

    plt.title("KM3NeT Angular Error")

    # plt.show()
    plt.savefig("KM3NeT Angular Error")

    return

plot_zenith_errors()
