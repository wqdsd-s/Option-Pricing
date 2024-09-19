import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

#qustion 2(a)

v0=0.1
alpha=0.45
beta=-5.105
gamma=0.3
S0=100
r=0.05
K=100

N=1000
T=5
'''
def q2a(rho):   
    step=1000
    delta=T/step

    #generate N stock paths
    Spath=np.zeros((N,step+1))
    Spath[:,0]=S0
    vpath=np.zeros((N,step+1))
    vpath[:,0]=v0
    for i in range(N):
        for j in range(1,step+1):
            dW=np.random.multivariate_normal([0,0],[[delta,rho*delta],[rho*delta,delta]])
            vpath[i,j]=vpath[i,j-1]+(alpha+beta*max(vpath[i,j-1],0))*delta+gamma*np.sqrt(max(vpath[i,j-1],0))*dW[1]
            Spath[i,j]=Spath[i,j-1]+r*Spath[i,j-1]*delta+Spath[i,j-1]*np.sqrt(max(vpath[i,j-1],0))*dW[0]

    #Lt and Ut series
    
    Lt=50*np.exp(0.138629*np.array(range(1001))*delta)
    Ut=200-50*np.exp(0.138629*np.array(range(1001))*delta)


    #check whether option is exercised for each path
    tauL=[]
    tauU=[]
    for i in range(N):
        try:
            tauL.append(np.argwhere(Spath[i]<=Lt)[0][0]*delta)

        except:
            tauL.append(T)

        try:
            tauU.append(np.argwhere(Spath[i]>=Ut)[0][0]*delta)
        except:
            tauU.append(T)
            
    check_tauL=np.array(tauL)<T    
    check_tauU=np.array(tauU)<T

    exercise=np.logical_or(check_tauL,check_tauU)

    N1=np.logical_and(exercise,np.array(tauL)<np.array(tauU))
    N1=len(N1[N1==True])
    N2=len(exercise[exercise==True])

    return N1/N2

for rho in [-0.74,0,0.74]:
    print(q2a(rho))
'''

#qustion 2(b)
def q2b(rho):   
    step=1000
    delta=T/step

    #generate N stock paths
    Spath=np.zeros((N,step+1))
    Spath[:,0]=S0
    vpath=np.zeros((N,step+1))
    vpath[:,0]=v0
    for i in range(N):
        for j in range(1,step+1):
            dW=np.random.multivariate_normal([0,0],[[delta,rho*delta],[rho*delta,delta]])
            vpath[i,j]=vpath[i,j-1]+(alpha+beta*max(vpath[i,j-1],0))*delta+gamma*np.sqrt(max(vpath[i,j-1],0))*dW[1]
            Spath[i,j]=Spath[i,j-1]+r*Spath[i,j-1]*delta+Spath[i,j-1]*np.sqrt(max(vpath[i,j-1],0))*dW[0]

    #Lt and Ut series
    
    Lt=50*np.exp(0.138629*np.array(range(1001))*delta)
    Ut=200-50*np.exp(0.138629*np.array(range(1001))*delta)


    #check whether option is exercised for each path
    tauL=[]
    tauU=[]
    tauL_i=[]
    tauU_i=[]
    for i in range(N):
        try:
            tauL.append(np.argwhere(Spath[i]<=Lt)[0][0]*delta)
                                
        except:
            tauL.append(T)
                   
        try:
            tauU.append(np.argwhere(Spath[i]>=Ut)[0][0]*delta)
        except:
            tauU.append(T)

        try:
            tauL_i.append(np.argwhere(Spath[i]<=Lt)[0][0])
        except:
            tauL_i.append(1000)
        try:
            tauU_i.append(np.argwhere(Spath[i]>=Ut)[0][0])
        except:
            tauU_i.append(1000)
            

    tauL_i=np.array(tauL_i)
    tauU_i=np.array(tauU_i)

    
    check_tauL=np.array(tauL)<T    
    check_tauU=np.array(tauU)<T

    exercise=np.logical_or(check_tauL,check_tauU)

    which1=np.logical_and((np.array(tauL)<=np.array(tauU)),exercise)
    which2=np.logical_and((np.array(tauL)>np.array(tauU)),exercise)

    index1=np.array(range(1000))[which1]
    i1=tauL_i[which1]

    index2=np.array(range(1000))[which2]
    i2=tauU_i[which2]

    payoff1=np.sum((K-Spath[index1,i1])*np.exp(-r*i1*delta))
    payoff2=np.sum((Spath[index2,i2]-K)*np.exp(-r*i2*delta))

    return (payoff1+payoff2)/1000

for rho in [-0.74,0,0.74]:
    print(q2b(rho))

        
    







         
