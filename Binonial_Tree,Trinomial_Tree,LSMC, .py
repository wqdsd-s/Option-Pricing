import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm


#LECTURE-4
#1.(a)

t=0.5
r=0.055
sigma=0.25
S0=180
K=170

#Pricing function
def Ame_Put_a(t,r,sigma,S0,K,n):
    delta=t/n
    c=(np.exp(-r*delta)+np.exp((r+sigma**2)*delta))/2
    d=c-np.sqrt(c**2-1)
    u=1/d
    p=(np.exp(r*delta)-d)/(u-d)

    S_path=np.zeros((n+1,n+1))
    S_path[:]=np.nan

    V_path=np.zeros((n+1,n+1))
    V_path[:]=np.nan

    #simulate for S in each node
    for i in range(n+1):
        for j in range(i+1):
            S_path[j,i]=S0*(u**(i-j))*(d**j)

    #calculate terminal option value
    V_path[:,-1]=np.maximum(K-S_path[:,-1],0)

    #calculate option value at each node backwards
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            V_path[j,i]=max((V_path[j,i+1]*p+V_path[j+1,i+1]*(1-p))*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]


price=[]

for n in [20,40,80,100,200,500]:
    price.append(Ame_Put_a(t,r,sigma,S0,K,n))

    
#1.(b)
    
#Pricing function
def Ame_Put_b(t,r,sigma,S0,K,n):
    delta=t/n
    u=np.exp((r-sigma**2/2)*delta+sigma*np.sqrt(delta))
    d=np.exp((r-sigma**2/2)*delta-sigma*np.sqrt(delta))
    p=1/2

    S_path=np.zeros((n+1,n+1))
    S_path[:]=np.nan

    V_path=np.zeros((n+1,n+1))
    V_path[:]=np.nan

    #simulate for S in each node
    for i in range(n+1):
        for j in range(i+1):
            S_path[j,i]=S0*(u**(i-j))*(d**j)

    #calculate terminal option value
    V_path[:,-1]=np.maximum(K-S_path[:,-1],0)
    
    #calculate option value at each node backwards
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            V_path[j,i]=max((V_path[j,i+1]*p+V_path[j+1,i+1]*(1-p))*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

#plot the prices in one graph   
price1=[]    
for n in [20,40,80,100,200,500]:
    price1.append(Ame_Put_b(t,r,sigma,S0,K,n))

plt.plot([20,40,80,100,200,500],price,label='a')   
plt.plot([20,40,80,100,200,500],price1,label='b')      
plt.xlabel("n")
plt.ylabel("Price")
plt.legend()
plt.show()


#2.(i)
t=0.5
r=0.055
sigma=0.25
K=170
n=500

#CRR Model function
def Ame_Put_CRR(t,r,sigma,S0,K,n):
    delta=t/n
    u=np.exp(sigma*np.sqrt(delta))
    d=1/u
    p=(np.exp(r*delta)-d)/(u-d)

    S_path=np.zeros((n+1,n+1))
    S_path[:]=np.nan

    V_path=np.zeros((n+1,n+1))
    V_path[:]=np.nan

    for i in range(n+1):
        for j in range(i+1):
            S_path[j,i]=S0*(u**(i-j))*(d**j)

    V_path[:,-1]=np.maximum(K-S_path[:,-1],0)

    for i in range(n-1,-1,-1):
        for j in range(i+1):
            V_path[j,i]=max((V_path[j,i+1]*p+V_path[j+1,i+1]*(1-p))*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

#Delta function
def Delta(t,r,sigma,S0,K,n,mu=0.15):
    return (Ame_Put_CRR(t,r,sigma,S0+mu,K,n)-Ame_Put_CRR(t,r,sigma,S0,K,n))/mu

#plot delta
delta_list1=[]
for S0 in range(170,192,2):
    delta_list1.append(Delta(t,r,sigma,S0,K,n,mu=0.15))
    
plt.plot(range(170,192,2),delta_list1)      
plt.xlabel("S0")
plt.ylabel("delta")
plt.show()


#2.(ii)
t=0.5
r=0.055
sigma=0.25
S0=180
K=170
n=500

delta_list2=[]
for t in np.arange(0,0.183,0.003):
    delta_list2.append(Delta(t,r,sigma,S0,K,n,mu=0.15))
    
plt.plot(np.arange(0,0.183,0.003),delta_list2)      
plt.xlabel("t")
plt.ylabel("delta")
plt.show()


#2.(iii)
t=0.5
r=0.055
sigma=0.25
S0=180
K=170
n=500

def Theta(t,r,sigma,S0,K,n,mu=0.15):
    return -(Ame_Put_CRR(t+mu,r,sigma,S0,K,n)-Ame_Put_CRR(t,r,sigma,S0,K,n))/mu

theta_list=[]
for t in np.arange(0,0.183,0.003):
    theta_list.append(Theta(t,r,sigma,S0,K,n,mu=0.15))
    
plt.plot(np.arange(0,0.183,0.003),theta_list)      
plt.xlabel("t")
plt.ylabel("theta")
plt.show()


#2.(iv)
t=0.5
r=0.055
sigma=0.25
S0=180
K=170
n=500

def Vega(t,r,sigma,S0,K,n,mu=0.15):
    return (Ame_Put_CRR(t,r,sigma+mu,S0,K,n)-Ame_Put_CRR(t,r,sigma,S0,K,n))/mu

vega_list=[]
for S0 in range(170,192,2):
    vega_list.append(Vega(t,r,sigma,S0,K,n,mu=0.15))
    
plt.plot(range(170,192,2),vega_list)      
plt.xlabel("S0")
plt.ylabel("vega")
plt.show()


#3.(a)
t=0.5
r=0.055
sigma=0.25
S0=180
K=170

#Trinomial Tree Pring Function
def Ame_Put_Tri_a(t,r,sigma,S0,K,n):
    delta=t/n
    
    d=np.exp(-sigma*np.sqrt(3*delta))
    u=1/d
    
    pd=(r*delta*(1-u)+(r*delta)**2+sigma**2*delta)/(u-d)/(1-d)
    pu=(r*delta*(1-d)+(r*delta)**2+sigma**2*delta)/(u-d)/(u-1)
    pm=1-pu-pd
    
    S_path=np.zeros((2*n+1,n+1))
    S_path[:]=np.nan

    V_path=np.zeros((2*n+1,n+1))
    V_path[:]=np.nan
    
    #simulate for S in each node
    for i in range(n+1):
        for j in range(2*i+1):
            if j<=i:
                S_path[j,i]=S0*(u**(i-j))
            else:
                S_path[j,i]=S0*(d**(j-i))

    #calculate terminal option value
    V_path[:,-1]=np.maximum(K-S_path[:,-1],0)

    #calculate option value at each node backwards
    for i in range(n-1,-1,-1):
        for j in range(2*i+1):
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

price3=[]    
for n in [20, 40, 70, 80, 100, 200, 500]:
    price3.append(Ame_Put_Tri_a(t,r,sigma,S0,K,n))

plt.plot([20, 40, 70, 80, 100, 200, 500],price3)      
plt.xlabel("n")
plt.ylabel("Price")
plt.show()


#3.(b)
t=0.5
r=0.055
sigma=0.25
S0=180
K=170

#Trinomial Tree Pring Function
def Ame_Put_Tri_b(t,r,sigma,S0,K,n):
    X0=np.log(S0)

    delta=t/n
    Xu=sigma*np.sqrt(3*delta)
    temp=r-sigma**2/2
    
    pd=((sigma**2*delta+temp**2*(delta**2))/(Xu**2)-temp*delta/Xu)/2
    pu=((sigma**2*delta+temp**2*(delta**2))/(Xu**2)+temp*delta/Xu)/2
    pm=1-pu-pd
    
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

price4=[]    
for n in [20, 40, 70, 80, 100, 200, 500]:
    price4.append(Ame_Put_Tri_b(t,r,sigma,S0,K,n))

plt.plot([20, 40, 70, 80, 100, 200, 500],price4)      
plt.xlabel("n")
plt.ylabel("Price")
plt.show()



#LECTURE-5
#4.(a)
def LSMC_Laguerre(k,T):
    S0=180
    sigma=0.25
    r=0.055
    K=170
    N=100000
    step=int(np.sqrt(N)*T)
    delta=T/step

    #generate N stock paths
    rand=(r-sigma**2/2)*delta+sigma*np.sqrt(delta)*np.random.normal(size=(N,step))
    log_incre=np.cumsum(rand,1)
    path=np.concatenate([(np.ones(N)*S0).reshape(-1,1),S0*np.exp(log_incre)],1)

    index=np.zeros((N,step+1))

    index[np.argwhere(path[:,-1]<K).flatten(),-1]=1

    for i in range(step-1,0,-1):
        atm_index=np.argwhere(path[:,i]<K).flatten()
        X=path[atm_index,i]

        #calculate for Y
        path_mask=path[atm_index,i+1:]
        value_mask=np.maximum(K-path_mask,0)
        index_mask=index[atm_index,i+1:]
        discount_mask=np.ones((len(atm_index),step-i))*np.exp(np.arange(1,step-i+1)*(-r)*delta)

        Y=(value_mask*index_mask*discount_mask).sum(1).reshape(-1,1)

        #calculate for L(x) with different k
        if k==2:
            LX=np.array([np.exp(-X/2),np.exp(-X/2)*(1-X)]).T
        elif k==3:
            LX=np.array([np.exp(-X/2),np.exp(-X/2)*(1-X),np.exp(-X/2)*(1-2*X+X**2/2)]).T
        elif k==4:
            LX=np.array([np.exp(-X/2),np.exp(-X/2)*(1-X),np.exp(-X/2)*(1-2*X+X**2/2),np.exp(-X/2)*(1-3*X+(X**2)*3/2-X**3/6)]).T
        elif k==5:
            LX=np.array([np.exp(-X/2),np.exp(-X/2)*(1-X),np.exp(-X/2)*(1-2*X+X**2/2),np.exp(-X/2)*(1-3*X+(X**2)*3/2-X**3/6),np.exp(-X/2)*(1-4*X+(X**2)*3-X**3/3*2+X**4/24)]).T
        
        #estimate for ECV
        try:
            a_est=np.linalg.inv((LX.T).dot(LX)).dot((LX.T).dot(Y))
            ECV=LX.dot(a_est).flatten()
            
        except:
            ECV=Y.flatten()
            
        #compare for EV and ECV and reset index path
        EV=K-X
        index[atm_index[EV>ECV],i]=1
        index[atm_index[EV>ECV],i+1:]=0

    #calculate final value at t=0
    path_mask=path[:,1:]
    value_mask=np.maximum(K-path_mask,0)
    index_mask=index[:,1:]
    discount_mask=np.ones((N,step))*np.exp(np.arange(1,step+1)*(-r)*delta)
    V=(value_mask*index_mask*discount_mask).sum()/N
    return V


#get all the 8 prices
res=pd.DataFrame()
for k in [2,3,4,5]:
    li=[]
    for T in [0.5,1.5]:
        li.append(LSMC_Laguerre(k,T))
    res[k]=li
res.columns.name='k'
res.index=[0.5,1.5]
res.index.name='T'
print(res)



#4.(b)
def LSMC_Hermite(k,T):
    S0=180
    sigma=0.25
    r=0.055
    K=170
    N=100000
    step=int(np.sqrt(N)*T)
    delta=T/step

    #generate N stock paths
    rand=(r-sigma**2/2)*delta+sigma*np.sqrt(delta)*np.random.normal(size=(N,step))
    log_incre=np.cumsum(rand,1)
    path=np.concatenate([(np.ones(N)*S0).reshape(-1,1),S0*np.exp(log_incre)],1)

    index=np.zeros((N,step+1))

    index[np.argwhere(path[:,-1]<K).flatten(),-1]=1

    for i in range(step-1,0,-1):
        atm_index=np.argwhere(path[:,i]<K).flatten()
        X=path[atm_index,i]

        #calculate for Y
        path_mask=path[atm_index,i+1:]
        value_mask=np.maximum(K-path_mask,0)
        index_mask=index[atm_index,i+1:]
        discount_mask=np.ones((len(atm_index),step-i))*np.exp(np.arange(1,step-i+1)*(-r)*delta)

        Y=(value_mask*index_mask*discount_mask).sum(1).reshape(-1,1)

        #calculate for L(x) with different k
        if k==2:
            LX=np.array([np.ones(len(X)),2*X]).T
        elif k==3:
            LX=np.array([np.ones(len(X)),2*X,4*X**2-2]).T
        elif k==4:
            LX=np.array([np.ones(len(X)),2*X,4*X**2-2,8*X**3-12*X]).T
        elif k==5:
            LX=np.array([np.ones(len(X)),2*X,4*X**2-2,8*X**3-12*X,16*X**4-56*X**2+16]).T

        #estimate for ECV
        try:
            a_est=np.linalg.inv((LX.T).dot(LX)).dot((LX.T).dot(Y))
            ECV=LX.dot(a_est).flatten()
            
        except:
            ECV=Y.flatten()
            
        #compare for EV and ECV and reset index path
        EV=K-X
        index[atm_index[EV>ECV],i]=1
        index[atm_index[EV>ECV],i+1:]=0

    #calculate final value at t=0
    path_mask=path[:,1:]
    value_mask=np.maximum(K-path_mask,0)
    index_mask=index[:,1:]
    discount_mask=np.ones((N,step))*np.exp(np.arange(1,step+1)*(-r)*delta)
    V=(value_mask*index_mask*discount_mask).sum()/N
    return V


#get all the 8 prices
res=pd.DataFrame()
for k in [2,3,4,5]:
    li=[]
    for T in [0.5,1.5]:
        li.append(LSMC_Hermite(k,T))
    res[k]=li
res.columns.name='k'
res.index=[0.5,1.5]
res.index.name='T'
print(res)


#4.(c)
def LSMC_Simple(k,T):
    S0=180
    sigma=0.25
    r=0.055
    K=170
    N=100000
    step=int(np.sqrt(N)*T)
    delta=T/step

    #generate N stock paths
    rand=(r-sigma**2/2)*delta+sigma*np.sqrt(delta)*np.random.normal(size=(N,step))
    log_incre=np.cumsum(rand,1)
    path=np.concatenate([(np.ones(N)*S0).reshape(-1,1),S0*np.exp(log_incre)],1)

    index=np.zeros((N,step+1))

    index[np.argwhere(path[:,-1]<K).flatten(),-1]=1

    for i in range(step-1,0,-1):
        atm_index=np.argwhere(path[:,i]<K).flatten()
        X=path[atm_index,i]

        #calculate for Y
        path_mask=path[atm_index,i+1:]
        value_mask=np.maximum(K-path_mask,0)
        index_mask=index[atm_index,i+1:]
        discount_mask=np.ones((len(atm_index),step-i))*np.exp(np.arange(1,step-i+1)*(-r)*delta)

        Y=(value_mask*index_mask*discount_mask).sum(1).reshape(-1,1)

        #calculate for L(x) with different k
        if k==2:
            LX=np.array([np.ones(len(X)),X]).T
        elif k==3:
            LX=np.array([np.ones(len(X)),X,X**2]).T
        elif k==4:
            LX=np.array([np.ones(len(X)),X,X**2,X**3]).T
        elif k==5:
            LX=np.array([np.ones(len(X)),X,X**2,X**3,X**4]).T

        #estimate for ECV
        try:
            a_est=np.linalg.inv((LX.T).dot(LX)).dot((LX.T).dot(Y))
            ECV=LX.dot(a_est).flatten()
            
        except:
            ECV=Y.flatten()
            
        #compare for EV and ECV and reset index path
        EV=K-X
        index[atm_index[EV>ECV],i]=1
        index[atm_index[EV>ECV],i+1:]=0

    #calculate final value at t=0
    path_mask=path[:,1:]
    value_mask=np.maximum(K-path_mask,0)
    index_mask=index[:,1:]
    discount_mask=np.ones((N,step))*np.exp(np.arange(1,step+1)*(-r)*delta)
    V=(value_mask*index_mask*discount_mask).sum()/N
    return V


#get all the 8 prices
res=pd.DataFrame()
for k in [2,3,4,5]:
    li=[]
    for T in [0.5,1.5]:
        li.append(LSMC_Simple(k,T))
    res[k]=li
res.columns.name='k'
res.index=[0.5,1.5]
res.index.name='T'
print(res)


#LECTURE-6 
#5.(a)
t=0.5
r=0.055
sigma=0.25
K=170
delta=0.002

def EFD(t,r,sigma,S0,K,delta,delta_X_multiplier):
    X0=np.log(S0)

    n=int(t/delta)
    Xu=sigma*np.sqrt(delta_X_multiplier*delta)
    temp=r-sigma**2/2
    
    pd=delta/2*(sigma**2/(Xu**2)-temp/Xu)
    pu=delta/2*(sigma**2/(Xu**2)+temp/Xu)
    pm=1-delta*(sigma**2)/(Xu**2)-r*delta
    
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]


#5.(b)
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]


#5.(c)
def CNFD(t,r,sigma,S0,K,delta,delta_X_multiplier):
    X0=np.log(S0)

    n=int(t/delta)
    Xu=sigma*np.sqrt(delta_X_multiplier*delta)
    temp=r-sigma**2/2
    
    pu=-delta/4*(sigma**2/(Xu**2)+temp/Xu)
    pd=-delta/4*(sigma**2/(Xu**2)-temp/Xu)
    pm=1+delta*(sigma**2)/(Xu**2)/2+r*delta/2
    
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

#summarize all methods
for delta_X_multiplier in [1,3,4]:
    list1=[]
    for S0 in range(170,191,1):
        list1.append(EFD(t,r,sigma,S0,K,delta,delta_X_multiplier))
    plt.plot(range(170,191,1),list1,label='EFD,dX=sigma*sqrt('+str(delta_X_multiplier)+'*dt)')
        
plt.plot(range(170,191,1),[Ame_Put_CRR(t,r,sigma,S0,K,250) for S0 in range(170,191,1)],label='bino_CRR')
plt.plot(range(170,191,1),[Ame_Put_Tri_a(t,r,sigma,S0,K,250) for S0 in range(170,191,1)],label='trino_3(a)')
plt.plot(range(170,191,1),[Ame_Put_Tri_b(t,r,sigma,S0,K,250) for S0 in range(170,191,1)],label='trino_3(b)')   
plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()


for delta_X_multiplier in [1,3,4]:
    list2=[]
    list3=[]
    for S0 in range(170,191,1):
        list2.append(IFD(t,r,sigma,S0,K,delta,delta_X_multiplier))
        list3.append(CNFD(t,r,sigma,S0,K,delta,delta_X_multiplier))
    plt.plot(range(170,191,1),list2,label='IFD,dX=sigma*sqrt('+str(delta_X_multiplier)+'*dt)')
    plt.plot(range(170,191,1),list3,label='CNFD,dX=sigma*sqrt('+str(delta_X_multiplier)+'*dt)')

plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()


#6.(a)
t=0.5
r=0.055
sigma=0.25
K=170
delta=0.002

def EFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier):
    X0=np.log(S0)

    n=int(t/delta)
    Xu=sigma*np.sqrt(delta_X_multiplier*delta)
    temp=r-sigma**2/2
    
    pd=delta/2*(sigma**2/(Xu**2)-temp/Xu)
    pu=delta/2*(sigma**2/(Xu**2)+temp/Xu)
    pm=1-delta*(sigma**2)/(Xu**2)-r*delta
    
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]


#5.(b)
def IFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier):
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]


#5.(c)
def CNFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier):
    X0=np.log(S0)

    n=int(t/delta)
    Xu=sigma*np.sqrt(delta_X_multiplier*delta)
    temp=r-sigma**2/2
    
    pu=-delta/4*(sigma**2/(Xu**2)+temp/Xu)
    pd=-delta/4*(sigma**2/(Xu**2)-temp/Xu)
    pm=1+delta*(sigma**2)/(Xu**2)/2+r*delta/2
    
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
            V_path[j,i]=max((V_path[j,i+1]*pu+V_path[j+1,i+1]*pm+V_path[j+2,i+1]*pd)*np.exp(-r*delta),K-S_path[j,i])

    return V_path[0,0]

#plot option price using three methods
delta_X_multiplier=3
list1=[]
list2=[]
list3=[]
for S0 in range(170,191,1):
    list1.append(EFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier))
    list2.append(IFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier))
    list3.append(CNFD_S(t,r,sigma,S0,K,delta,delta_X_multiplier))
    
plt.plot(range(170,191,1),list1,label='EFD')
plt.plot(range(170,191,1),list2,label='IFD')
plt.plot(range(170,191,1),list3,label='CNFD')
plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()

plt.plot(range(170,191,1),list1,label='EFD')     
plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()

plt.plot(range(170,191,1),list2,label='IFD')
plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()

plt.plot(range(170,191,1),list3,label='CNFD')
plt.xlabel("S0")
plt.ylabel("Price")
plt.legend()
plt.show()


