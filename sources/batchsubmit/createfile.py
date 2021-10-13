import pandas as pd

# first define some parameters
xlo = 0
xhi = 10
xstep = 1
ylo = 0
yhi = 10
ystep = 1

xlen = int((xhi - xlo + 1)/xstep)
ylen = int((yhi - ylo + 1)/ystep)

# create list of DM^2 for the dict
x = []
y = []

for xi in range(xlen):
    for yi in range(ylen):
        x.append(xlo + xi * xstep)
        y.append(ylo + yi * ystep)
        yi += 1
    xi += 1

d = {'X' : x, 'Y' : y}

df = pd.DataFrame(d)

print(df.loc[df['X'] == 1])
