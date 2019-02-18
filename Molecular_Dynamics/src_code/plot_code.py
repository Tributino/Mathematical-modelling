#plots for assignment 
#example of matplotlib usage

import scipy
import math
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

df=pd.read_csv("N100")

scipy.nan
fig, ax=plt.subplots()
plt.xlim(0.0001,0.005)
plt.ylim(0,0.8)
ax.scatter(1/df['Volume_1'], df['Psim'], color="r", label = "Simulated Pressure")
ax.plot(1/df['Volume_1'],df['Pthe'],label="Gas with Excluded Volume - Theoretical Pressure", color="orange")
ax.plot(1/df['Volume_1'], df['Pideal'], label="Ideal Gas - Theoretical Pressure", color="gray")
ax.set(xlabel="Volume^(-1) (1/V)", ylabel="Pressure (P)", title = "Particles Based Modeling: Pressure vs Volume")
ax.grid()
plt.legend()
fig.savefig("volume_inverso.png")
