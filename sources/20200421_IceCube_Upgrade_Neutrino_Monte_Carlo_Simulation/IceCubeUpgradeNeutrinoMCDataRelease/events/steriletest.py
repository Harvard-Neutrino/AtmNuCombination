import sterileprob
# import sterileopt
import numpy as np
import time

# from sterileopt import t2kNH, novaNH
from sterileprob import prob
import matplotlib.pyplot as plt

start = time.time()
for i in range(1000):
	prob(2, 1, "nova", "NH", "t2k", m4 = 1, t14 = 0.17, t24= 0.143486, t34= 0, \
	mand13 = True, vald13 = -1.539785, manE = False, valE = 1, d14 = -0.004956, anti = False, avg = True)
end = time.time()

print("used time", end - start)

# print(prob(2, 1, "t2k", "NH", m4 = 1, t14 = 0., t24= 0., t34= 0, \
# 	mand13 = True, vald13 = np.pi, d14 = 0, d34 = 0, anti = False))
# print(prob(2, 1, "nova", "NH", m4 = 1, t14 = 0., t24= 0., t34= 0, \
	# mand13 = False, vald13 = np.pi/4, d14 = 0, d34 = 0, anti = True))
# print(t2kNH(0.15, 0.15, -np.pi/4, np.pi/4))
# print(prob(2, 1, "t2k", "NH", m4 = 1, t14 = 0.1704502851872749, t24= 0.18367733438823333, t34= 0, \
# 	mand13 = True, vald13 = -0.06855384197330401, d14 = 1.5197067798275956, d34 = 0, anti = False))
# print(prob(2, 1, "t2k", "NH", m4 = 1, t14 = 0.1704502851872749, t24= 0.18367733438823333, t34= 0, \
# 	mand13 = True, vald13 = -0.06855384197330401, d14 = 1.5197067798275956, d34 = 0, anti = True))
# sterileopt.NE()

# print(prob(2, 1, "nova", "NH", m4 = 1, t14 = 0., t24= 0., t34= 0, \
# 	mand13 = False, vald13 = np.pi/4, d14 = 0, d34 = 0, anti = False))
# print(prob(2, 1, "nova", "NH", m4 = 1, t14 = 0., t24= 0., t34= 0, \
# 	mand13 = False, vald13 = np.pi/4, d14 = 0, d34 = 0, anti = True))
# print(prob(2, 1, "t2k", "NH", "t2k", m4 = 1, t14 = 0.15, t24= 0.1, t34= 0, \
#  	mand13 = True, vald13 = 0, d14 = 0, d34 = 0, anti = False))
# def PlotMatterAvg(savename = "blank2.png"):
#     x = np.linspace(0, 1, 10000)
#     y = np.zeros(len(x))
#     z = np.zeros(len(x))
#     for i in range(1, len(x)):
#         y[i] = prob(2, 1, "t2k", "NH", "t2k", m4 = 1, t14 = 0.15, t24= 0.1, t34= 0, \
# 						  	mand13 = True, vald13 = 0, manE = True, valE = x[i], d14 = 0, d34 = 0, anti = False, avg = False)
#         z[i] = prob(2, 1, "t2k", "NH", "t2k", m4 = 1, t14 = 0.15, t24= 0.1, t34= 0, \
# 						  	mand13 = True, vald13 = 0, manE = True, valE = x[i], d14 = 0, d34 = 0, anti = False, avg = True)
#     plt.plot(x,y)
#     plt.plot(x,y)
#     plt.grid()
#     plt.savefig(savename)

# print(prob(2, 1, "t2k", "NH", "t2k", m4 = 1, t14 = 0.1, t24= 0.1, t34 = 0.1, \
# 	mand13 = False, vald13 = 3.7, manE = False, valE = 0.2, d14 = 0, d34 = 0, anti = False, avg = False))
# print(prob(2, 1, "nova", "NH", "t2k", m4 = 1, t14 = 0.17, t24= 0.143486, t34= 0, \
# 	mand13 = True, vald13 = -1.539785, manE = False, valE = 1, d14 = -0.004956, anti = False, avg = True))

# print(prob(2, 1, "t2k", "NH", "t2k", m4 = 1, t14 = 0.15, t24= 0.1, t34= 0, \
# 	mand13 = False, vald13 = 0, manE = True, valE = 0.2, d14 = 0, d34 = 0, anti = False, avg = True))