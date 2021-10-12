#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:18:47 2021

@author: tercio
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import log as ln
from math import e



class Receptors:
    """
    tsd - time of synaptic delay
    tr  - time of EPSP/IPSP rise
    td - time of EPSP/IPSP decay
    """
    
    def __init__(self, Amax=0, tsd=0, tr=0, td=0):
        self.__Amax = Amax;
        self.__tsd = tsd
        self.__tr = tr
        self.__td = td
        
    
    def AMPA(self):
        
       self.__Amax = 5; 
       self.__tsd = 1
       self.__tr = 2
       self.__td = 13
               
       return [self.__Amax, self.__tsd, self.__tr, self.__td]
   
    def NMDA (self):
        
       self.__Amax = 1; 
       self.__tsd = 1
       self.__tr = 2
       self.__td = 13
               
       return [self.__Amax, self.__tsd, self.__tr, self.__td]
    def GABA (self):
        
       self.__Amax = -2.5; 
       self.__tsd = 1
       self.__tr = 2
       self.__td = 10
               
       return [self.__Amax, self.__tsd, self.__tr, self.__td]
    
 
def synaptic (t, Amax, tsd, tr, td ):
    """ Calculates the typical time courses of excitatory or inhibitory
    postsynaptic potentials
    
    tsd — time of synaptic delay
    tr — time of EPSP/IPSP rise
    td — time of EPSP/IPSP decay, 1 millisecond (ms) = 2 steps of i
    """
    
    if t == tsd:
        x = 0
        
    if t > tsd  and t <= tr:
        x = (Amax/tr)*(t-tsd)
        
    if t > tr and t <= td: 
        x = (Amax/td)*((td-(t-(tr+tsd))))
        
          
        
    return x

def adaptation (S):
    """  The adjust on typical time postsynaptic potentials
    S(i) is the actual value of summarised potential in the k compartment, 
    ReP = –80mV - resting potential value
    Rsp = –10 mV - reverse synaptic potential

    """
    ReP = -80;
    Rsp = -10;
    result = ((ReP-Rsp) - (ReP-S))/(ReP-Rsp)
    
    return result
     
def weight(LSW,k, NE):
    """ Weight of compartment
    LSW  — weight of the most distal dendrite input.
    NE = 13 — a number of excitatory inputs.
    """
    
        
    result = ((1-LSW)/NE-1)*(k-NE)+1
    return result

def function_influence (LSW, N, a, b):
    """Influence of a compartment on b compartment (for excitatory inputs)
    N = total number of inputs
    a, b є N
    """
    
    if b >= a:
        x = 1 - ((1-LSW)/(N-1))*(a-b)
    if (a-b) <= (N/4):
        x = 1 - (N/4)*(a-b)
    else:
        x = 0
    

def summarised_potential(E):
    """Summarised potential
    NE — number of excitatory inputs, 
    Inf (m,k) - Influence of a compartment on b compartment 
    E(m;0) the  actual value of the appropriate register.
    ReP - resting potential value
    """
    ReP = -80;     
    y = ReP + (E - ReP)#*function_influence(m,k))
    
    
    return y
    
def mem (C):
    """Long-term synaptic potentiation (LTP) 
    Memory (Mem)
    C(k; i) time of memory for compartment
    clog parameter = 2.3026."""
    clog = 2.3026
    x = 1 - ln ((C+1)/6*clog)
    
    return x
    
def power (M):
    """modell that phenomenological LTP event
    powerA = 9 is a parameter,
    M(k;i) — actual value of SF(i) for EPSP(NMDA) in appropriate register."""
    powerA = 5e-5;
    ReP = -80;
   
    return powerA*(M-ReP)

def C (C,power):
    """Time of memory duration
    """
    x = np.zeros(len(C))
    x[0] = C[0]
    for i in range(len(C)-1):
        x[i+1] = x[i] + e**(power[i]) - 1
    return x
                   


##### Paramestros iniciais ####

# “k” is the number of the dendritic compartment,
# “i” is a number of an area in a particular register table.


#Um compartilmento e uma sinápse,  k = 1, e i = 1
ReP = -80;
k = 1
NE = 1
NI = 1
LSW = 0.2
EPSPd = 4.5
PQ = 10

ampa = Receptors().AMPA()
nmda = Receptors().NMDA()
gaba = Receptors().GABA()

t = np.zeros(20)

m = np.zeros(20)

for i in range(len(t)):
    #t[i] = random.randint(1, 10)
    if i < 10:
        t[i]=3
    else:
        t[i]=3
        

for i in range(len(m)):
    if i < 10:
        m[i]=1
    else:
        m[i]=10

 
E = np.zeros(len(t))
E [0] = synaptic (t[0], ampa[0], ampa[1],ampa[2], ampa[3])

M = np.zeros(len(t))
M [0] = synaptic (t[0], nmda[0], nmda[1],nmda[2], nmda[3])

I = np.zeros(len(t))
I [0] = synaptic (t[0], gaba[0], gaba[1],gaba[2], gaba[3])
        
for i in range(len(t)-1):
    
    sf = synaptic (t[i], ampa[0], ampa[1],ampa[2], ampa[3])    
    E [i+1] = E [i]+adaptation(ReP+sf)*mem(m[i])*sf
    
    sf = synaptic (t[i], nmda[0], nmda[1],nmda[2], nmda[3])    
    M[i+1] = M [i]+sf
    
    sf = synaptic (t[i], gaba[0], gaba[1],gaba[2], gaba[3])    
    I[i+1] = I [i]+sf
    
    



PPS = np.zeros(len(t))

for i in range (len(t)):
    """Summarised postsynaptic potential in neuron"""

    PPS[i] = ReP + ((weight(LSW, k,NE)*(E[i]-ReP))+(I[i]-ReP))
        
 
Mpower = np.zeros(len(M))

for i in range(len(M)):
    Mpower[i] = power(M[i])
    
    
timeMemory = C(M, Mpower)    
        
plt.plot(E, 'k')
plt.plot(M, 'b')
plt.plot(I, 'r')

plt.plot(timeMemory)

x = np.arange(1,21,1)

plt.plot(x,PPS, 'g')
plt.ylabel("EPSP (mv)")
        
    




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


x = np.arange(1,21,1)
y = PPS



fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False)
#ax.set_aspect('equal')
line, = ax.plot([], [], 'or', lw=2)
ax.set_xlim(0, np.max(x))
ax.set_ylim(np.min(y), np.max(y)+10)
ax.set_ylabel("EPSP (mv)")
ax.set_xlabel("Trials")



def animate(i):
    thisx =  x[0:i]
    thisy = y[0:i]

    line.set_data(thisx, thisy)
    
    return line

ani = ani.FuncAnimation(fig, animate, np.arange(1, len(y)),)
ani.save('EPSP.mp4', fps=10)
plt.show()











