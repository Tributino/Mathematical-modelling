#!/usr/bin/env python
# coding: utf-8

# # Euler method
# 
# Euler's method is a numerical method to solve first order first degree differential equation with a given initial value. 
# 
# 
# \begin{equation*}
#     \begin{cases}
#     \cfrac{dy}{dt} = f(t,y),
#     \\
#     y(0) = y_0.
#     \end{cases}
# \end{equation*}
# 
# It is the most basic explicit method for numerical integration of ordinary differential equations.
# 
# The idea is that while the curve is initially unknown, its starting point, $P_0=(t_0,y_0)$ is known. Then, from the differential equation, the slope to the curve at $P_0$ can be computed, and so, the tangent line.
# 
# Take a small step along that tangent line up to a point $P_1 = (t_1,y_1)$. Along this small step, the slope does not change too much, so $P_1$ will be close to the curve. After several steps, a polygonal curve $P_{0}, P_{1}, P_{2}, P_{3}...$ is computed. In general, this curve does not diverge too far from the original unknown curve, and the error between the two curves can be made small if the step size is small enough and the interval of computation is finite.
# 
# <img src="graph_euler_method.png" style="width:400px">
# 
# Than, from the fundamental theorem of calculus we can conclude that for the general differential equation
# 
# $$\frac{dy}{dt} = f(t,y),$$
# 
# the $n$-th step of Euler's method is given by
# 
# $$y((n+1)h) = y(nh) + hf(t,y(nh)),$$
# 
# where $h = t_{n+1}-t_n$ is some time step choosed.
# 

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Let's try the Euler Method in the below equation for differents uniform step size $h$.
# 
# \begin{equation*}
#     \begin{cases}
#     \cfrac{dy}{dt} = y - \frac{1}{2}e^{\frac{t}{2}}\sin(5t) + 5e^{\frac{t}{2}}\cos(5t),
#     \\
#     y(0) = 0.
#     \end{cases}
# \end{equation*}
# 
# For more examples see [Paul's Online Notes](http://tutorial.math.lamar.edu/Classes/DE/EulersMethod.aspx)

# In[3]:


def euler(step):
    
    # Initializations
    h = step                                # step 
    y_0 = 0                                 # initial value
    t_0 = 0                                 # initial time
    t_end = 5                               # stopping time
    n_steps = int(round((t_end-t_0)/h))     # total number of steps
    
    # Let's create a dictionary with the poits P_n=(y_n,t_n), where t_n's are keys and y_n's are values
    
    P_dic = {}                              
    P_dic[t_0] = y_0                          # add the first point P_0=(y_0,t_0)

    # Euler's method
    t = t_0
    y = P_dic[t_0]

    for i in range (1, n_steps + 1):
        dydt = y-1/2*np.exp(t/2)*np.sin(5*t) + 5*np.exp(t/2)*np.cos(5*t)        # calculate the derivative dy/dt 
        P_dic[i*h] = y + h*dydt
        t = i*h                                                                 # Next Point 
        y = P_dic[i*h] 

    return P_dic


# In[4]:


for step in [0.5, 0.05, 0.005]:

    t = np.copy(list(euler(step).keys()))
    y = np.copy(list(euler(step).values()))
   
    # Plot the results
    
    plt.plot(t, y, linewidth = 2, linestyle='dashed')       # plot P_n's


y_exact = np.exp(t/2)*np.sin(5*t)       # Analytical Solution
plt.plot(t, y_exact, linewidth = 2)     # plot y_exact vs. t
   
plt.plot([0,6],[0,0],'k',linewidth=1)
plt.legend(('for h = 0.5','for h = 0.05', 'for h = 0.005', 'y_exact','y = 0'))    
plt.xlabel('time', fontsize = 20)
plt.ylabel('y(t)', fontsize = 20)

plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.grid(True)                          # show grid 
plt.axis([0, 6, -15, 25])               # define the axes
plt.savefig("contrast.png")             # saving image
plt.show()                              # show the plot


# ### Let's now create a Dataframe with time as  index and columns that compare the analytic solution with the numerical solution.

# In[5]:


# To create a array with the y aproximate by Euler Method 
def aprox(step, times):
    aprox_arr = []
    for t in times:
        y = euler(step)[t]
        aprox_arr.append(y)
    return aprox_arr

# To create a array with the exact y 
def exact(times):
    times = np.array(times)
    exact = np.exp(times/2)*np.sin(5*times)
    return exact

# To know how close the approximation and the exact solution are? 
def Per_Errors(step, times):
    error = 100*((np.abs(exact(times) - aprox(step, times)))/np.abs(exact(times)))
    return error


# In[11]:


# Creating the DataFrame.
times = [1, 2, 3, 4, 5]
d = {'Time': times, 'Exact': exact(times), 
     'h = 0.5': aprox(0.5, times), 'h = 0.05': aprox(0.05, times),'h = 0.005': aprox(0.005, times),
     "Error_for h = 0.5": Per_Errors(0.5, times),
     "Error_for h = 0.05": Per_Errors(0.05, times),
     "Error_for h = 0.005": Per_Errors(0.005, times)}
df = pd.DataFrame(data=d)
df


# In[12]:


#Set the DataFrame index using the columns Time.
df = df.set_index('Time')


# In[13]:


# Round a DataFrame to a variable number of decimal places.
df = df.round(2)
df


# In[14]:


# To Apply the simble %
df['Error_for h = 0.5'] = df['Error_for h = 0.5'].astype(str) + '%'
df['Error_for h = 0.05'] = df['Error_for h = 0.05'].astype(str) + '%'
df['Error_for h = 0.005'] = df['Error_for h = 0.005'].astype(str) + '%'
df


# In[75]:


# Saving the Dataframe as a csv
df.to_csv('euler_method_ex.csv')


# # References
# 
# - [Mechatronical Modelling/Chapter 7](http://moodle.autolab.uni-pannon.hu/Mecha_tananyag/mechatronikai_modellezes_angol/ch07.html#d0e9679)
# - [Paul's Online Notes](http://tutorial.math.lamar.edu/Classes/DE/EulersMethod.aspx)
# - [Course Mathematical Modelling Basics](https://www.edx.org/course/mathematical-modelling-basics)
# - [Wikipedia/Euler method](https://en.wikipedia.org/wiki/Euler_method)
#    
# 
