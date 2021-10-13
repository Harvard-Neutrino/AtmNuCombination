import sys
import pandas as pd

def foo(x, y):
    return x ** 2 + y ** 2

pkl = sys.argv[1]
i = int(sys.argv[2])

unpkl = pd.read_pickle(pkl)

unpklX = unpkl['X']
unpklY = unpkl['Y']

currentX = unpklX.iloc[i]
currentY = unpklY.iloc[i]

print(foo(currentX, currentY))
