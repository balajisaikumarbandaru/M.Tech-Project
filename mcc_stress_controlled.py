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
def mcc_stresscontrolled(t,p_0v,p_eq_0,v_0,M,l,k,G,mu,T,p0_0,dqdt):
    t,p_0v,p_eq_0,v_0,M,l,k,G,mu,T,p0_0
    dt=float(input("enter dt"))
    days=int(input("enter no of days"))
    t=int(days/dt)
    p=np.zeros((t,1))
    q=np.zeros((t,1))
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
    d_v=np.zeros((t,1))
    v=np.zeros((t,1))
    dp_0=np.zeros((t,1))
    p_0=np.zeros((t,1))
    eta=np.zeros((t,1))
    d_eta=np.zeros((t,1))
    t1=np.zeros((t,1))
    dq=dqdt*dt
    dp=dq/3
    p[0]=p_0v;p_eq[0]=p_eq_0;p_0[0]=p0_0;v[0]=1.73977
    for i in range(1,t):
            p[i]=p[i-1]+dp
            q[i]=q[i-1]+dq
            eta[i]=q[i]/p[i]
            d_eta[i]=eta[i]-eta[i-1]
            p_eq[i]=((M**2+(eta[i])**2)/(M**2))*p[i]
            d_epsilonve[i]=(k/v[i-1])*(dp/p[i])
            d_epsilonse[i]=dq/(3*G)
            d_epsilonvpdt[i]=((M**2-(eta[i])**2)/(M**2))*(mu/(v[i-1]*T))*(p_eq[i]/p_0[i-1])**((l-k)/mu)
            d_epsilonvp[i]=d_epsilonvpdt[i]*dt
            d_epsilonsvpdt[i]=((2*(eta[i])**2)/(M**2))*(mu/(v[i-1]*T))*(p_eq[i]/p_0[i-1])**((l-k)/mu)
            d_epsilonsvp[i]= d_epsilonsvpdt[i]*dt
            d_epsilons[i]=d_epsilonsvp[i]+d_epsilonse[i]
            d_epsilonv[i]=d_epsilonvp[i]+d_epsilonve[i]
            epsilonv[i]=epsilonv[i-1]+d_epsilonv[i]
            epsilons[i]=epsilons[i-1]+d_epsilons[i]
            d_v[i]=(v[i-1])*epsilonv[i]
            v[i]=v[i-1]-d_v[i]
            dp_0[i]=(p_0[i-1]*d_epsilonv[i]*(v[i]))/(l-k)
            p_0[i]=p_0[i-1]+dp_0[i]
            if q[i]>p[i]:
                break
    return t1,p_0,dp_0,d_v,p,q,v,d_v,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons
t1,p_0,dp_0,d_v,p,q,v,d_v,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons=mcc_stresscontrolled(60,200,200,1.73977,0.9,0.15,0.03,10000,0.006,1,300,3)
p=p[p!=0]
q=removezeros(p,q)
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
d_v=removezeros(p,d_v)
v=removezeros(p,v)
dp_0=removezeros(p,dp_0)
p_0=removezeros(p,p_0)
t1=removezeros(p,t1)
fig=plt.figure(1)
plt.plot(p,v)

t1,p_0,dp_0,d_v,p,q,v,d_v,p_eq,d_epsilonve,d_epsilonse,d_epsilonvpdt,d_epsilonvp,epsilons,d_epsilonsvp,d_epsilonv,epsilonv,d_epsilonsvpdt,d_epsilons=mcc_stresscontrolled(60,200,200,1.73977,0.9,0.15,0.03,10000,0.006,1,300,30)
p=p[p!=0]
q=removezeros(p,q)
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
d_v=removezeros(p,d_v)
v=removezeros(p,v)
dp_0=removezeros(p,dp_0)
p_0=removezeros(p,p_0)
t1=removezeros(p,t1)
plt.plot(p,v)