import os
try:
    import numpy as np
    from matplotlib import pyplot as plt
except:
    os.system('cmd/k "pip install numpy"')
    os.system('cmd/k "pip install matplotlib"')
    import numpy as np
    from matplotlib import pyplot as plt
"""
funtions
"""
def MCC(p_0,e_0,M,l,k,dp):
    pf=(3*M*p_0)/(3-M)
    #iter=int((pf-p_0)/dp)
    iter=100
    dp=(pf-p_0)/iter
    p=np.zeros((iter,1))
    q=np.zeros((iter,1))
    e=np.zeros((iter,1))
    eta=np.zeros((iter,1))
    d_eta=np.zeros((iter,1))
    d_epsilonv=np.zeros((iter,1))
    epsilonv=np.zeros((iter,1))
    d_e=np.zeros((iter,1))
    e=np.zeros((iter,1))
    d_epsilons=np.zeros((iter,1))
    epsilons=np.zeros((iter,1))
    epsilona=np.zeros((iter,1))
    d_v=np.zeros((iter,1))
    v=np.zeros((iter,1))
    e[0]=e_0;p[1]=p_0;
    eta[1]=q[1]/p[1]
    d_eta[1]=eta[1]-eta[0]
    d_epsilonv[1]=(l/(1+e[1]))*((dp/p[1])+(1-k/l)*(2*eta[1]*d_eta[1])/(M**2+eta[1]**2))
    d_e[1]=(1+e[0])*d_epsilonv[1]
    e[1]=e[0]-d_e[1]
    v[0]=1.73977
    for i in range(2,iter):
            p[i]=p[i-1]+dp
            q[i]=q[i-1]+3*dp
            eta[i]=q[i]/p[i]
            d_eta[i]=eta[i]-eta[i-1]
            d_epsilonv[i]=(l/(1+e[i]))*((dp/p[i])+(1-k/l)*(2*eta[i]*d_eta[i])/(M**2+eta[i]**2))
            d_e[i]=(1+e[i-1])*d_epsilonv[i]
            d_epsilons[i]=((2*eta[i])/(M**2-eta[i]**2))*((l-k)/(1+e[i-1]))*((dp)/p[i]+(2*eta[i]*d_eta[i])/(M**2+eta[i]**2))
            epsilons[i]=epsilons[i-1]+d_epsilons[i]
            epsilonv[i]=epsilonv[i-1]+d_epsilonv[i]
            epsilona[i]= epsilonv[i]/3+epsilons[i]
            e[i]=e[i-1]-d_e[i]
            d_v[i]=(v[i-1])*epsilonv[i]
            v[i]=v[i-1]-d_v[i]
            if q[i]>(3*M*p[1])/(3-M):
                break
                print(i)
    return epsilons,epsilona,q,e,epsilonv,p,v
"""
data set 01 - p_o=[200,150,100], e_0 = 0.889, M=1,lambda=0.174,k=0.026
data set 02 - p_o=[300], e_0 = 0.889, M=0.9,lambda=0.15,k=0.03
"""
plt.figure(1)
epsilons,epsilona,q,e,epsilonv,p,v=MCC(200,.889,1,.174,.026,7.2)
dict1={"epsilons":epsilons,"epsilona":epsilona,"q":q,"e":e,"epsilonv":epsilonv,"p":p}
plt.plot(epsilona,q)

epsilons,epsilona,q,e,epsilonv,p,v=MCC(150.8,.889,1,.174,.026,7.2)
plt.plot(epsilona,q)
dict2={"epsilons":epsilons,"epsilona":epsilona,"q":q,"e":e,"epsilonv":epsilonv,"p":p}

epsilons,epsilona,q,e,epsilonv,p,v=MCC(100,.889,1,.174,.026,7.2)
plt.plot(epsilona,q)
dict3={"epsilons":epsilons,"epsilona":epsilona,"q":q,"e":e,"epsilonv":epsilonv,"p":p}
plt.ylabel("Devitoric stress")
plt.xlabel("axial strain")
plt.legend(["p_o=200kPa","p_o=150kPa","p_o=100kPa"])
plt.figure(2)
epsilons,epsilona,q,e,epsilonv,p,v=MCC(200,.889,1,.174,.026,7.2)
plt.plot(epsilona,epsilonv)
epsilons,epsilona,q,e,epsilonv,p,v=MCC(150,.889,1,.174,.026,7.2)
plt.plot(epsilona,epsilonv)
epsilons,epsilona,q,e,epsilonv,p,v=MCC(100,.889,1,.174,.026,7.2)
plt.plot(epsilona,epsilonv)
plt.ylabel("Volumetric strain")
plt.xlabel("axial strain")
plt.legend(["p_o=200kPa","p_o=150kPa","p_o=100kPa"])
