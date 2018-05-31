import math
import numpy as np
import plotly
import plotly.graph_objs as go
from scipy.integrate import odeint

file = open('C:/Users/Polina/Desktop/универ/vm/ex1.txt', 'r') 
lines = file.readlines() 
ex1 = np.zeros((21,81))
for i in range(21):
    ex1[i] = list(map(lambda x: float(x),lines[i].split(',')))
ex2 = np.zeros((21,21))

#start conditions
T = 0.1
def a(w): return 1
def b(w): return 0
def c(point): return math.sin(point[0])
def alpha1(t):return 1
def alpha2(t):return 1
def beta1(t):return 1
def beta2(t):return 1

#depends on u(x,t)
#w = (x,t)
def u(w): #return w[0] + w[1]
    return w[0]**3 + w[1]**3 
def dudx(w): #return 1
    return 3*w[0]**2 
def dudt(w): #return 1
    return 3*w[1]**2 
def du2dx(w): #return 0
    return 6*w[0] 

def fi(x): return u([x,0])
def f(w): return dudt(w) - du2dx(w) - u(w)*math.sin(w[0])
def alpha(t): return alpha1(t)*u([0,t]) - alpha2(t)*dudx([0,t]) 
def beta(t): return beta1(t)*u([1,t]) + beta2(t)*dudx([1,t])

def rnd(arr,n):
    for i in range(len(arr)):
        for j in range(len(arr[0])):
            arr[i][j] = round(arr[i][j],n)
    return(arr)

def sweep(n,A,B,C,G):
    S,T = [],[]
    S.append(C[0]/B[0])
    T.append(-G[0]/B[0])
    for i in range(1,n+1):
        S.append(C[i]/(B[i]-A[i]*S[i-1]))
        T.append((A[i]*T[i-1]-G[i])/(B[i]-A[i]*S[i-1]))
    Y = [0 for i in range(n+1)]
    Y[-1] = T[-1]
    for i in range(n-1,-1,-1):
        Y[i] = S[i]*Y[i+1]+T[i]
    return Y

def setka(N,M):
    W = np.zeros((N+1,M+1,2))
    for i in range(N+1):
        for k in range(M+1):
            W[i][k] = (i/N,k*T/M)
    return W

def Lh(i,k,W,U):
    h = 1/(len(W)-1)
    return a(W[i][k])*(U[i+1][k]-2*U[i][k]+U[i-1][k])/h**2 + b(W[i][k])*(U[i+1][k]-U[i-1][k])/(2*h) + c(W[i][k])*U[i][k]

def solution(N,M,W):
    U = np.zeros((N+1,M+1))
    h, tau = 1/N, T/M
    for i in range(N+1): U[i][0] = fi(W[i][0][0])
    for k in range(1,M+1):
        for i in range(1,N): U[i][k] = U[i][k-1] + tau*(Lh(i,k-1,W,U)+f(W[i][k-1]))
        U[0][k] = (alpha(W[0][k][1]) + alpha2(W[0][k][1])*(4*U[1][k]-U[2][k])/(2*h)) / (alpha1(W[0][k][1])+3*alpha2(W[0][k][1])/(2*h))
        U[N][k] = (beta(W[0][k][1]) - beta2(W[0][k][1])*(-4*U[N-1][k]+U[N-2][k])/(2*h)) / (beta1(W[0][k][1])+3*beta2(W[0][k][1])/(2*h))
    return U

def get_M(N,T):
    h = 1/N
    tau = h**2/2
    return int(T/tau)+1

def weights_schema(N,M,sigma):
    W = setka(N,M)
    U = np.zeros((N+1,M+1))
    h, tau = 1/N, T/M
    for i in range(N+1): U[i][0] = fi(W[i][0][0])
    if sigma == 0:
        for k in range(1,M+1):         
            for i in range(1,N): U[i][k] = U[i][k-1] + tau*(Lh(i,k-1,W,U)+f(W[i][k-1]))
            U[0][k] = (alpha(W[0][k][1])+alpha2(W[0][k][1])*U[1][k]/h) / (alpha1(W[0][k][1])+alpha2(W[0][k][1])/h)
            U[N][k] = (beta(W[0][k][1])+beta2(W[0][k][1])*U[N-1][k]/h) / (beta1(W[0][k][1])+beta2(W[0][k][1])/h)
    elif sigma == 1/2:
        for k in range(1,M+1):
            A,B,C,G = [],[],[],[]
            A.append(0)
            B.append(-alpha1(W[0][k][1])-alpha2(W[0][k][1])/h)
            C.append(-alpha2(W[0][k][1])/h)
            G.append(alpha(W[0][k][1]))
            for i in range(0,N):
                A.append(sigma*a(W[i][k])/(h**2) - b(W[i][k])/(2*h))
                B.append(2*sigma*a(W[i][k])/(h**2) - c(W[i][k])+ 1/tau)
                C.append(sigma*a(W[i][k])/(h**2) + b(W[i][k])/(2*h))
                G.append(-U[i][k-1]/tau - (1-sigma)*Lh(i,k-1,W,U) - f([W[i][k][0],W[i][k][1]-tau/2]))    
            A.append(-beta2(W[0][k][1])/h)
            B.append(-beta1(W[0][k][1])-beta2(W[0][k][1])/h)
            C.append(0)
            G.append(beta(W[0][k][1]))   
            for i in range(0,N+1): U[i][k] = sweep(N,A,B,C,G)[i]
    return U 



def table1(N,M,W,U):
    W2 = setka(5,5)
    U2 = np.zeros((6,6))
    for i in range(len(W)):
        for j in range(len(W[0])):
            if i%(N/5)==0 and j%(M/5)==0: U2[i//(N//5)][j//(M//5)] = U[i][j]
    v = np.concatenate((np.array([[W2[i][0][0] for i in range(len(W2))]]).T,rnd(U2,3)),axis=1).T
    trace = go.Table(
        header=dict(values=['x/t']+[str(W2[0][i][1]) for i in range(len(W2[0]))]),
        cells=dict(values= v.tolist()))    
    data = [trace] 
    plotly.offline.plot(data)    

def norm(u,smth):
    mx = 0
    k1 = (len(smth)-1)//(len(u)-1)
    k2 = (len(smth[0])-1)//(len(u[0])-1)
    for i in range(len(u)):
        for j in range(len(u[0])):
            if abs(u[i][j] - smth[i*k1][j*k2]) > mx:
                mx = abs(u[i][j] - smth[i*k1][j*k2])
    return mx

def table2(U_lst,U2_lst,num):
    if num == 1:
        norm1,norm2 = [norm(U_lst[i],ex1) for i in range(3)],[norm(U_lst[i],U2_lst[i]) for i in range(3)]
    if num == 2:
        norm1,norm2 = [0,0,0],[norm(U_lst[i],U2_lst[i]) for i in range(3)]
    trace = go.Table(
        header=dict(values=['h','tau','||ex - u(h,tau)||','||u(h,tau) - u(2h,tau2)||']),
        cells=dict(values=[[0.2,0.1,0.05],[0.02,0.005,0.00125],norm1,norm2]))    
    data = [trace] 
    plotly.offline.plot(data)    

def plot_table(N_lst,M,schema,table_num):
    U_lst = []
    U2_lst = []
    for i in range(len(N_lst)):
        N = N_lst[i]       
        if schema == 1: M = get_M(N,T)
        h,tau = 1/N,T/M
        W = setka(N,M)
        if schema == 1: U = solution(N,M,W)
        if schema == 2: U = weights_schema(N,M,1/2) #change sigma here
        U_lst.append(U)
        if i > 0: U2_lst.append(U)
        if i+1 == table_num: table1(N,M,W,U)
    N = N_lst[-1]*2       
    M = get_M(N,T)
    h,tau = 1/N,T/M
    W = setka(N,M)
    if schema == 1: U = solution(N,M,W)
    if schema == 2: U = weights_schema(N,M,1/2) #change sigma here
    U2_lst.append(U)
    if table_num == 0: table2(U_lst,U2_lst,1)

#plot_table([5,10,20],0,1,1) #schema1, N = 5
#plot_table([5,10,20],0,1,2) #schema1, N = 10
#plot_table([5,10,20],0,1,3) #schema1, N = 20
#plot_table([5,10,20],0,1,0) #schema1, table for all
#plot_table([5,10,20],10,2,1) #schema2, N = 5, M = 10
#plot_table([5,10,20],10,2,2) #schema2, N = 10, M = 10
#plot_table([5,10,20],10,2,3) #schema2, N = 20, M = 10
plot_table([5,10,20],10,2,0) #schema2, table for all
#plot_table([5,10,20],100,2,1) #schema2, N = 5, M = 100
#plot_table([5,10,20],100,2,2) #schema2, N = 10, M = 100
#plot_table([5,10,20],100,2,3) #schema2, N = 20, M = 100
#plot_table([5,10,20],100,2,0) #schema2, table for all

