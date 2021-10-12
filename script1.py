#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 18:47:24 2021

@author: tercio
"""

import matplotlib.pyplot as plt
import numpy as np
from random import randint
  
def AMPA (t):
    """ Calculates the typical time courses of excitatory or inhibitory
    postsynaptic potentials
    
    tsd — time of synaptic delay
    tr — time of EPSP/IPSP rise
    td — time of EPSP/IPSP decay, 1 millisecond (ms) = 2 steps of i
    """
    Amax, tsd, tr, td =  1, 1, 2, 13          
        
    return 0*(t==tsd) or ((Amax/tr)*(t-tsd))*(t > tsd and t <= tr) or ((Amax/td)*((td-(t-(tr+tsd)))))*(t > tr and t <= td)

def NMDA (t):
    """ Calculates the typical time courses of excitatory or inhibitory
    postsynaptic potentials
    
    tsd — time of synaptic delay
    tr — time of EPSP/IPSP rise
    td — time of EPSP/IPSP decay, 1 millisecond (ms) = 2 steps of i
    """
    Amax, tsd, tr, td =  1, 1, 2, 13          
        
    return 0*(t==tsd) or ((Amax/tr)*(t-tsd))*(t > tsd and t <= tr) or ((Amax/td)*((td-(t-(tr+tsd)))))*(t > tr and t <= td)


def GABA (t):
    """ Calculates the typical time courses of excitatory or inhibitory
    postsynaptic potentials
    
    tsd — time of synaptic delay
    tr — time of EPSP/IPSP rise
    td — time of EPSP/IPSP decay, 1 millisecond (ms) = 2 steps of i
    """
    Amax, tsd, tr, td =  -2.5, 1, 2, 10          
        
    return 0*(t==tsd) or ((Amax/tr)*(t-tsd))*(t > tsd and t <= tr) or ((Amax/td)*((td-(t-(tr+tsd)))))*(t > tr and t <= td)



simulationTime = 100; #in milliseconds
deltaT=.01;
t=np.arange(0,simulationTime, deltaT);


"""specify the external current I==="""
changeTimes = [0]; #in milliseconds
I = np.zeros(len(t))

I[0:500] = [randint(0, 15) for p in range(0, 500)]; 
I[501:2000] = 0; 
I[2001:len(t)] = [randint(0, 15) for p in range(0, 7999)]
#Set externally applied current across time
#Here, first 500 timesteps are at current of 50, next 1500 timesteps at
#current of zero (resets resting potential of neuron), and the rest of
#timesteps are at constant current

ampa = np.zeros(len(t))
nmda = np.zeros(len(t))
gaba = np.zeros(len(t))

V= np.zeros(len(t))
b= np.zeros(len(t))
A = np.zeros(len(t))
S = np.zeros(len(t))
S1 = np.zeros(len(t))


ReP = -80; 
Rsp = -10
C = 1
deltaT=.01;
M = -.1

for i in range(len(t)-1):
    ampa[i] = AMPA(I[i])
    nmda[i] = NMDA(I[i])
    gaba[i] = GABA(I[i])
    soma = ampa[i]+nmda[i]+gaba[i]
    
    S [i] = ReP + (ampa[i] - ReP)
   
    
            
    A = ((ReP-Rsp) - (ReP-S[i]))/(ReP-Rsp)
    b[i+1] =  b[i]+A*M*AMPA(I[i])
    
    V[i+1] = V[i] + deltaT*A*soma/C;
    
     
   
fig, axs = plt.subplots(4)
axs[0].plot(t, ampa,'r', label='AMPA')
axs[1].plot(t, nmda,'g', label='NMDA')
axs[2].plot(t, gaba,'b', label='GABA')
axs[3].plot(t, I,'k', label='Applied current across time')

axs[0].legend()
axs[1].legend()
axs[2].legend()
axs[3].legend()

axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xticks([])

axs[0].set_ylabel('Voltage (mv)')
axs[1].set_ylabel('Voltage (mv)')
axs[2].set_ylabel('Voltage (mv)')
axs[3].set_ylabel('Voltage (mv)')

axs[3].set_xlabel('time (ms)')


plt.plot(t, b,'b')



