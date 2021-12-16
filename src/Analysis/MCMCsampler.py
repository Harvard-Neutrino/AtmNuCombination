import numpy as np
import emcee
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import corner

def osci(theta, phi, data):
    return data*np.sin(theta)*np.sin(theta)*np.cos(phi)*np.cos(phi)

def DiffChi(theta, dataraw):
    th, phi = theta
    O = osci(th,phi,dataraw)
    E  = osci(0.5,1.,dataraw)
    chi = (O-E)*(O-E)/E
    return -0.5*chi


data = np.random.randint(100,size=20)+10

ndim, nwalkers = 2, 20
p0 = np.random.rand(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, DiffChi, args=[data])
sampler.run_mcmc(p0, 1000)



samples = sampler.get_chain(flat=True)
th = samples - 2*3.14159265*(samples//(2*3.14159265))

fig = corner.corner(th);

# fig = plt.figure()
# gs = GridSpec(1, 1, wspace=0.06, hspace=0.1)
# a1 = plt.subplot(gs[0])
# th = samples[:, 0]
# th = th - 2*3.14159265*(th//(2*3.14159265))
# print(th.size)
# a1.hist(th, 100, color="k", histtype="step")
# a1.set_xlabel(r"$\theta_1$")
# a1.set_ylabel(r"$p(\theta_1)$")
plt.show()

