import math
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sterileutil

# Returns the all sterile-active mixings including the CP phase
def PMNS(t12, t13, t23, t14, t24, t34, d13, d14, d34):
	#print("All4MixCP Matrix generating...")
	U12 = np.array([[np.cos(t12), np.sin(t12), 0, 0],
					[-np.sin(t12), np.cos(t12), 0, 0],
					[0, 0, 1, 0],
					[0, 0, 0, 1]], dtype = complex)
	U14 = np.array([[np.cos(t14), 0, 0, np.sin(t14) * np.exp(-1j * d14)],
					[0, 1, 0, 0],
					[0, 0, 1, 0],
					[-np.sin(t14) * np.exp(1j * d14), 0, 0, np.cos(t14)]])
	U23 = np.array([[1, 0, 0, 0],
					[0, np.cos(t23), np.sin(t23), 0],
					[0, -np.sin(t23), np.cos(t23), 0],
					[0, 0, 0, 1]])
	U13 = np.array([[np.cos(t13), 0, np.sin(t13) * np.exp(-1j * d13), 0],
					[0, 1, 0 ,0],
					[-np.sin(t13) * np.exp(1j * d13), 0, np.cos(t13), 0],
					[0, 0, 0, 1]])
	U24 = np.array([[1, 0, 0, 0],
					[0, np.cos(t24), 0, np.sin(t24)],
					[0, 0, 1, 0],
					[0, -np.sin(t24), 0, np.cos(t24)]])
	U34 = np.array([[1, 0, 0, 0],
					[0, 1, 0, 0],
					[0, 0, np.cos(t34), np.sin(t34) * np.exp(-1j * d34)],
					[0, 0, -np.sin(t34) * np.exp(1j * d34), np.cos(t34)]])
	R1 = np.matmul(U34, U24)
	R2 = np.matmul(R1, U14)
	R3 = np.matmul(R2, U23)
	R4 = np.matmul(R3, U13)
	res = np.matmul(R4, U12)
	# print("generated U matrix")
	# print(res)
	return res
