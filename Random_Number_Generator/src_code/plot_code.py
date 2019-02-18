#plots for assignment 

import scipy
import math
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

N = 1000;
hist_bins=150;

df=pd.read_csv("outo.dat")
arg_dist = scipy.linspace(0,math.pi,100000)

#graphics specifications
#frequency at classes were normalized to fit at sin(x) distribution
#bins were positioned between 0 and pi
#ax plot is the distribution sin(x) between [0,pi]

fig, ax=plt.subplots()
plt.xlim(0,3.5)
plt.ylim(0,1.2)
ax.scatter(df['i']*math.pi/hist_bins, df['hist_b']/max(df['hist_b']), color="blue", label = "Simulated Numbers")
ax.plot(arg_dist,scipy.sin(arg_dist),label="Sin (x) Distribution", color="magenta")
ax.set(xlabel="x", ylabel=r"Prob of pick x $\rho$(x)", title = "Distribution Sinusoidal")
ax.grid()
plt.legend()
fig.savefig("assig_2.png")


