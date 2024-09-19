import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

#LECTURE-7
#1.(a)
V0=20000
L0=22000
miu=-0.1
sigma=0.2
gamma=-0.4
lambda1=0.2
T=5
r0=0.055
delta=0.25
lambda2=0.4
alpha=0.7
episilon=0.95

N=1000

def Proj3_1a(lambda1,T):
    deltat=T/N

    R=r0+delta*lambda2
    r=R/12
    n=T*12
    PMT=L0*r/(1-(1/((1+r)**n)))
    a=PMT/r
    b=PMT/r/((1+r)**n)
    c=1+r
    beta=(episilon-alpha)/T

    #qt and Lt series
    qt=alpha+beta*np.array(range(1001))*deltat
    Lt=a-b*(c**(12*np.array(range(1001))*deltat))

    #do 1000 simulation of value path
    Vpath=np.zeros((N,1001))
    Vpath[:,0]=V0
    
    for i in range(N):
        for j in range(1,1001):
            Vpath[i,j]=Vpath[i,j-1]+Vpath[i,j-1]*(miu*deltat+sigma*np.random.normal(0,deltat)+gamma*np.random.poisson(lambda1*deltat))

    #calculate stopping time tau and payoff for each path
    tau=[]
    payoff=[]
    for i in range(N):
        try:
            #get the stopping time tau of ith path
            tau_i=np.argwhere(Vpath[i]<=(qt*Lt))[0][0]
            tau.append(tau_i*deltat)

            #calculate discounted payoff at time tau
            payoff.append(max(Lt[tau_i]-episilon*Vpath[i,tau_i],0)*np.exp(-r0*tau_i*deltat))
        except:
            pass
    return np.mean(payoff)

print(Proj3_1a(0.2,5))

#plot for lambda1 changes
price1=[]
for l in np.linspace(0.05,0.4,8):
    price1.append(Proj3_1a(l,5))

plt.plot(np.linspace(0.05,0.4,8),price1)
plt.xlabel("lambda1")
plt.ylabel("Price")
plt.show()



#1.(b)
V0=20000
L0=22000
miu=-0.1
sigma=0.2
gamma=-0.4
lambda1=0.2
T=5
r0=0.055
delta=0.25
lambda2=0.4
alpha=0.7
episilon=0.95

N=1000

def Proj3_1b(lambda1,T):
    deltat=T/N

    R=r0+delta*lambda2
    r=R/12
    n=T*12
    PMT=L0*r/(1-(1/((1+r)**n)))
    a=PMT/r
    b=PMT/r/((1+r)**n)
    c=1+r
    beta=(episilon-alpha)/T
    
    qt=alpha+beta*np.array(range(1001))*deltat
    Lt=a-b*(c**(12*np.array(range(1001))*deltat))

    #do 1000 simulation of value path
    Vpath=np.zeros((N,1001))
    Vpath[:,0]=V0
    
    for i in range(N):
        for j in range(1,1001):
            Vpath[i,j]=Vpath[i,j-1]+Vpath[i,j-1]*(miu*deltat+sigma*np.random.normal(0,deltat)+gamma*np.random.poisson(lambda1*deltat))

    #calculate stopping time tau for each path
    tau=[]
    for i in range(N):
        try:
            tau.append(np.argwhere(Vpath[i]<=(qt*Lt))[0][0]*deltat)
        except:
            pass
    return len(tau)/N

print(Proj3_1b(0.2,5))

#plot for lambda1 changes
dp=[]
for l in np.linspace(0.05,0.4,8):
    dp.append(Proj3_1b(l,5))

plt.plot(np.linspace(0.05,0.4,8),dp)
plt.xlabel("lambda1")
plt.ylabel("default probability")
plt.show()



#1.(c)
V0=20000
L0=22000
miu=-0.1
sigma=0.2
gamma=-0.4
lambda1=0.2
T=5
r0=0.055
delta=0.25
lambda2=0.4
alpha=0.7
episilon=0.95

N=1000

def Proj3_1c(lambda1,T):
    deltat=T/N

    R=r0+delta*lambda2
    r=R/12
    n=T*12
    PMT=L0*r/(1-(1/((1+r)**n)))
    a=PMT/r
    b=PMT/r/((1+r)**n)
    c=1+r
    beta=(episilon-alpha)/T
    
    qt=alpha+beta*np.array(range(1001))*deltat
    Lt=a-b*(c**(12*np.array(range(1001))*deltat))

    #do 1000 simulation of value path
    Vpath=np.zeros((N,1001))
    Vpath[:,0]=V0
    
    for i in range(N):
        for j in range(1,1001):
            Vpath[i,j]=Vpath[i,j-1]+Vpath[i,j-1]*(miu*deltat+sigma*np.random.normal(0,deltat)+gamma*np.random.poisson(lambda1*deltat))

    #calculate stopping time tau for each path
    tau=[]
    for i in range(N):
        try:
            tau.append(np.argwhere(Vpath[i]<=(qt*Lt))[0][0]*deltat)
        except:
            pass
    return np.mean(tau)

print(Proj3_1c(0.2,5))

#plot for lambda1 changes
tau1=[]
for l in np.linspace(0.05,0.4,8):
    tau1.append(Proj3_1c(l,5))

plt.plot(np.linspace(0.05,0.4,8),tau1)
plt.xlabel("lambda1")
plt.ylabel("tau")
plt.show()



#2.(a)

v0=0.1
alpha=0.45
beta=-5.105
gamma=0.25
S0=100
r=0.05
rho=-0.75
K=100
T=1

N=1000
#do 1000 simulation of stock price path
Spath=np.zeros((N,1001))
Spath[:,0]=S0
vpath=np.zeros((N,1001))
vpath[:,0]=v0
for i in range(N):
    for j in range(1,1001):
        dW=np.random.multivariate_normal([0,0],[[0.001,rho*0.001],[rho*0.001,0.001]])
        vpath[i,j]=vpath[i,j-1]+(alpha+beta*max(vpath[i,j-1],0))*0.001+gamma*np.sqrt(max(vpath[i,j-1],0))*dW[1]
        Spath[i,j]=Spath[i,j-1]+r*Spath[i,j-1]*0.001+Spath[i,j-1]*np.sqrt(max(vpath[i,j-1],0))*dW[0]
        
#calculate option price
check=Spath<94
check=np.sum(check,1).astype(bool)
VT=np.maximum(K-Spath[:,-1],0)
VT[check]=0

print(np.exp(-r*T)*np.mean(VT))

#2.(b)
barrier=np.array(range(1001))*0.001*6+91
check=Spath<barrier
check=np.sum(check,1).astype(bool)
VT=np.maximum(K-Spath[:,-1],0)
VT[check]=0

print(np.exp(-r*T)*np.mean(VT))

#2.(c)
barrier=-np.array(range(1001))*0.001*6+97
check=Spath<barrier
check=np.sum(check,1).astype(bool)
VT=np.maximum(K-Spath[:,-1],0)
VT[check]=0

print(np.exp(-r*T)*np.mean(VT))



#LECTURE-8
#3.(a)

r0=0.05
sigma=0.12
k=0.92
rbar=0.055
T=4

N=1000
#do 1000 simulation of interest rate path
path=np.zeros((N,4001))
path[:,0]=r0
for i in range(N):
    for j in range(1,4001):
        path[i,j]=path[i,j-1]+k*(rbar-path[i,j-1])/1000+sigma*np.sqrt(path[i,j-1])*np.random.normal(0,0.001)

#calculate bond price        
price=0        
for t in range(1,9):
    mean=np.mean(np.exp(-np.sum(path[:,:t*500],1)*0.001))
    if t<8:
        price=price+30*mean
    else:
        price=price+1030*mean
        
print(price)



#3.(b)
t=0
T=0.5
S=1
K=980

#do 1000 simulation of interest rate path
path=np.zeros((N,1001))
path[:,0]=r0
for i in range(N):
    for j in range(1,1001):
        path[i,j]=path[i,j-1]+k*(rbar-path[i,j-1])/1000+sigma*np.sqrt(path[i,j-1])*np.random.normal(0,0.001)

#calculate call option price
bond_price=1000*np.exp(-0.001*np.sum(path[:,500:-1],1))
call_price=np.mean(np.exp(-0.001*np.sum(path[:,0:500],1))*np.maximum(bond_price-K,0))

print(call_price)


#3.(c)

def IFD(t,r,sigma,S0,K,delta,delta_X_multiplier):
    X0=np.log(S0)

    n=int(t/delta)
    Xu=sigma*np.sqrt(delta_X_multiplier*delta)
    temp=r-sigma**2/2
    
    pu=-delta/2*(sigma**2/(Xu**2)+temp/Xu)
    pd=-delta/2*(sigma**2/(Xu**2)-temp/Xu)
    pm=1+delta*(sigma**2)/(Xu**2)+r*delta
    
    S_path=np.zeros((2*n+1,n+1))
    S_path[:]=np.nan
    
    X_path=np.zeros((2*n+1,n+1))
    X_path[:]=np.nan

    V_path=np.zeros((2*n+1,n+1))
    V_path[:]=np.nan

    for i in range(n+1):
        for j in range(2*i+1):
                X_path[j,i]=X0+(i-j)*Xu

    S_path=np.exp(X_path)

    V_path[:,-1]=np.maximum(K-S_path[:,-1],0)

    for i in range(n-1,-1,-1):
        for j in range(2*i+1):
            V_path[j,i]=(V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta)

    return V_path[0,0]



#4
x0=0
y0=0
r0=0
rho=0.7
a=0.1
b=0.3
sigma=0.05
eta=0.09
c=0.055

K=950
T=0.5
S=1
N=1000

def putprice(rho):
    #do 1000 simulation of interest rate path rt
    path_x=np.zeros((N,1001))
    path_x[:,0]=x0
    path_y=np.zeros((N,1001))
    path_y[:,0]=y0
    path_r=np.zeros((N,1001))
    path_r[:,0]=r0
    for i in range(N):
        for j in range(1,1001):
            dW=np.random.multivariate_normal([0,0],[[0.001,rho*0.001],[rho*0.001,0.001]])
            path_x[i,j]=path_x[i,j-1]-a*path_x[i,j-1]*0.001+sigma*dW[0]
            path_y[i,j]=path_y[i,j-1]-b*path_y[i,j-1]*0.001+eta*dW[1]
            path_r[i,j]=path_x[i,j]+path_y[i,j]+c

    #calculate put option price
    bond_price=1000*np.exp(-0.001*np.sum(path_r[:,500:-1],1))
    put_price=np.mean(np.exp(-0.001*np.sum(path_r[:,0:500],1))*np.maximum(K-bond_price,0))

    return put_price

print(putprice(0.7))

#plot for rho changes
price4=[]
for rho in np.linspace(-0.7,0.7,15):
    price4.append(putprice(rho))

plt.plot(np.linspace(-0.7,0.7,15),price4)
plt.xlabel("rho")
plt.ylabel("Price")
plt.show()


#5.(a)
WAC=0.08
V=100000
r0=0.078
k=0.6
rbar=0.08
sigma=0.12
T=30

N=1000

def MBS(k,rbar,sigma):
    #do 1000 simulation of interest rate path
    path=np.zeros((N,3001))
    path[:,0]=r0
    for i in range(N):
        for j in range(1,3001):
            path[i,j]=path[i,j-1]+k*(rbar-path[i,j-1])*0.01+sigma*np.sqrt(path[i,j-1])*np.random.normal(0,0.01)
            
    cf=np.array([[V*WAC]*29+[V*WAC+V]]*N)
    df=np.exp(-np.cumsum(path,1)*0.01)[:,99:-1:100]
    return np.mean(np.sum(cf*df,1))
    
print(MBS(k,rbar,sigma))

#plot for rbar changes
price5=[]
for rb in np.linspace(0.04,0.1,7):
    price5.append(MBS(k,rb,sigma))

plt.plot(np.linspace(0.04,0.1,7),price5)
plt.xlabel("rbar")
plt.ylabel("Price")
plt.show()

#5.(b)
Pbar=98000
def MBS_const(r0):
    r1=np.array([r0]*30)
    cf=np.array([V*WAC]*29+[V*WAC+V])
    df=np.exp(-np.cumsum(r1))
    return np.sum(cf*df)

#use binary method to calculate const interest rate and OAS
def OAS(k,rbar,sigma,Pbar):
    a=0.06
    b=0.09
    for i in range(100):
        c=(a+b)/2
        if MBS_const(c)>=Pbar:
            a=c
        else:
            b=c
    return c-r0
    
print(OAS(k,rbar,sigma,Pbar))


#5.(c)
Pbar=98000
def MBS_const(r0):
    r1=np.array([r0]*30)
    cf=np.array([V*WAC]*29+[V*WAC+V])
    df=np.exp(-np.cumsum(r1))
    return np.sum(cf*df)

def OAS(k,rbar,sigma,Pbar):
    a=0.06
    b=0.09
    for i in range(100):
        c=(a+b)/2
        if MBS_const(c)>=Pbar:
            a=c
        else:
            b=c
    return c-r0
    
print(OAS(k,rbar,sigma,Pbar))




