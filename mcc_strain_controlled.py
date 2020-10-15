# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 21:16:23 2020

@author: Balaji Sai Kumar B
"""
import os
try:
    import numpy as np
    from matplotlib import pyplot as plt
except:
    os.system('cmd/k "pip install numpy"')
    os.system('cmd/k "pip install matplotlib"')
    import numpy as np
    from matplotlib import pyplot as plt
def removezeros(a,b):
    lent=len(a)
    b=b[0:lent]
    return b
def mcc_strain_controlled(p_0v,p_eq_0,v_0,M,l,k,G,mu,T,p0_0):
    dt=float(input("enter dt")) #select dt as 0.0001
    desdt=float(input("enter strain rate"))
    t=int(2/dt)
    p=np.zeros((t,1))
    d_p=np.zeros((t,1))
    q=np.zeros((t,1))
    d_q=np.zeros((t,1))
    p_eq=np.zeros((t,1))
    d_epsilonve=np.zeros((t,1))
    d_epsilonse=np.zeros((t,1))
    d_epsilonvpdt=np.zeros((t,1))
    d_epsilonvp=np.zeros((t,1))
    d_epsilonsvp=np.zeros((t,1))
    d_epsilonsvpdt=np.zeros((t,1))
    d_epsilons=np.zeros((t,1))
    d_epsilonv=np.zeros((t,1))
    epsilonv=np.zeros((t,1))
    epsilons=np.zeros((t,1))
    epsilona=np.zeros((t,1))
    d_v=np.zeros((t,1))
    v=np.zeros((t,1))
    dp_0=np.zeros((t,1))
    p_0=np.zeros((t,1))
    t1=np.zeros((t,1))
    eta=np.zeros((t,1))
    d_eta=np.zeros((t,1))
    p[0] = p_0v;v[0]=v_0;p_eq[0]=p_eq_0;p_0[0]=p0_0;eta[0]=q[0]/p[0]
    for i in range(1,t):
        t1[i]=t1[i-1]+dt
        epsilons[i]=epsilons[i-1]+desdt*dt
        d_epsilons[i]=epsilons[i]-epsilons[i-1]
        d_epsilonsvpdt[i]=((2*(eta[i-1])**2)/(M**2))*(mu/(v[i-1]*T))*(p_eq[i-1]/p_0[i-1])**((l-k)/mu)
        d_epsilonsvp[i]= d_epsilonsvpdt[i]*dt
        d_epsilonse[i]=d_epsilons[i]-d_epsilonsvp[i]
        d_q[i]=d_epsilonse[i]*3*G
        q[i]=q[i-1]+d_q[i]
        p[i]=p[i-1]+d_q[i]/3
        eta[i]=q[i]/p[i]
        d_p[i]=p[i]-p[i-1]
        d_epsilonve[i]=(k/v[i-1])*(d_p[i]/p[i])
        d_epsilonvpdt[i]=((M**2-(eta[i])**2)/(M**2))*(mu/(v[i-1]*T))*(p_eq[i]/p_0[i-1])**((l-k)/mu)
        d_epsilonv[i]=d_epsilonvp[i]+d_epsilonve[i]
        epsilonv[i]=epsilonv[i-1]+d_epsilonv[i]
        d_v[i]=(v[i-1])*epsilonv[i]
        v[i]=v[i-1]-d_v[i]
        dp_0[i]=(p_0[i-1]*d_epsilonv[i]*(v[i]))/(l-k)
        p_0[i]=p_0[i-1]+dp_0[i]
        epsilona[i]= epsilonv[i]/3+epsilons[i]
        if p[i]<q[i]:
            break
    return epsilona,t1,p_0,dp_0,d_v,p,q,v,d_v,d_q,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons
epsilona,t1,p_0,dp_0,d_v,p,q,v,d_v,d_q,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons=mcc_strain_controlled(200,200,1.73977,0.9,0.15,0.03,10000,0.006,1,300)
p=p[p!=0]
q=removezeros(p,q)
d_q=removezeros(p,d_q)
p_eq=removezeros(p,p_eq)
d_epsilonve=removezeros(p,d_epsilonve)
d_epsilonse=removezeros(p,d_epsilonse)
d_epsilonvpdt=removezeros(p,d_epsilonvpdt)
d_epsilonvp=removezeros(p,d_epsilonvp)
d_epsilonsvp=removezeros(p,d_epsilonsvp)
d_epsilonsvpdt=removezeros(p,d_epsilonsvpdt)
d_epsilons=removezeros(p,d_epsilons)
d_epsilonv=removezeros(p,d_epsilonv)
epsilonv=removezeros(p,epsilonv)
epsilons=removezeros(p,epsilons)
epsilona=removezeros(p,epsilona)
d_v=removezeros(p,d_v)
v=removezeros(p,v)
dp_0=removezeros(p,dp_0)
p_0=removezeros(p,p_0)
t1=removezeros(p,t1)
fig=plt.figure(1)
plt.plot(epsilona,q)

epsilona,t1,p_0,dp_0,d_v,p,q,v,d_v,d_q,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons=mcc_strain_controlled(200,200,1.73977,0.9,0.15,0.03,10000,0.006,1,300)
p=p[p!=0]
q=removezeros(p,q)
d_q=removezeros(p,d_q)
p_eq=removezeros(p,p_eq)
d_epsilonve=removezeros(p,d_epsilonve)
d_epsilonse=removezeros(p,d_epsilonse)
d_epsilonvpdt=removezeros(p,d_epsilonvpdt)
d_epsilonvp=removezeros(p,d_epsilonvp)
d_epsilonsvp=removezeros(p,d_epsilonsvp)
d_epsilonsvpdt=removezeros(p,d_epsilonsvpdt)
d_epsilons=removezeros(p,d_epsilons)
d_epsilonv=removezeros(p,d_epsilonv)
epsilonv=removezeros(p,epsilonv)
epsilons=removezeros(p,epsilons)
epsilona=removezeros(p,epsilona)
d_v=removezeros(p,d_v)
v=removezeros(p,v)
dp_0=removezeros(p,dp_0)
p_0=removezeros(p,p_0)
t1=removezeros(p,t1)
plt.plot(epsilona,q)