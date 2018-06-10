import math
import numpy as np
import plotly
import plotly.graph_objs as go
from scipy.integrate import odeint

#start conditions
T = 0.1
def a(x,t): return 1
def b(x,t): return 0
def c(x,t): return math.sin(x)
def alpha1(t):return 1
def alpha2(t):return 1
def beta1(t):return 1
def beta2(t):return 1

#depends on u(x,t)
def u(x,t):
    #return x+t
    return x**3 + t**3 
def dudx(x,t):
    #return 1
    return 3*x**2 
def dudt(x,t):
    #return 1
    return 3*t**2 
def du2dx(x,t):
    #return 0
    return 6*x 

def fi(x): return u(x,0)
def f(x,t): return dudt(x,t) - du2dx(x,t) - u(x,t)*math.sin(x)
def alpha(t): return alpha1(t)*u(0,t) - alpha2(t)*dudx(0,t) 
def beta(t): return beta1(t)*u(1,t) + beta2(t)*dudx(1,t) 

def rnd(arr,n): #rounds all values in 2d arr
    return([[round(arr[i][j],n) for j in range(len(arr[0]))] for i in range(len(arr))])

def sweep(n,A,B,C,G): #system range = n+1
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

def grid(N,M): #N,M - numbers of parts
    W = np.zeros((N+1,M+1,2))
    for i in range(N+1):
        for k in range(M+1):
            W[i][k] = (i/N,k*T/M)
    return W

#values of exact solution in points 6x6
ex_sol = np.zeros((6,6))
W_tmp = grid(5,5)
ex_sol = [[u(W_tmp[i][j][0],W_tmp[i][j][1]) for j in range(6)] for i in range(6)]

def Lh(i,k,h,tau,U):
    return a(i*h,tau*k)*(U[i+1][k]-2*U[i][k]+U[i-1][k])/(h**2) + b(i*h,tau*k)*(U[i+1][k]-U[i-1][k])/(2*h) + c(i*h,tau*k)*U[i][k]

def explicit_scheme(N,M):
    U = np.zeros((N+1,M+1))
    h,tau = 1/N,T/M
    for i in range(N+1): U[i][0] = fi(i*h)
    for k in range(1,M+1):
        for i in range(1,N): U[i][k] = U[i][k-1] + tau*(Lh(i,k-1,h,tau,U)+f(h*i,tau*(k-1)))
        U[0][k] = (alpha(tau*k) + alpha2(tau*k)*(4*U[1][k]-U[2][k])/(2*h)) / (alpha1(tau*k)+3*alpha2(tau*k)/(2*h))
        U[N][k] = (beta(tau*k) - beta2(tau*k)*(-4*U[N-1][k]+U[N-2][k])/(2*h)) / (beta1(tau*k)+3*beta2(tau*k)/(2*h))
    return U

def get_M(N,T): #calculate M for explicit scheme
    h = 1/N
    tau = h**2/2
    return int(T/tau)+1

def weights_scheme(N,M,sigma):
    U = np.zeros((N+1,M+1))
    h,tau = 1/N,T/M
    for i in range(N+1): U[i][0] = fi(i*h)
    if sigma == 0:
        for k in range(1,M+1):         
            for i in range(1,N): U[i][k] = U[i][k-1] + tau*(Lh(i,k-1,h,tau,U) + f(i*h,tau*(k-1)))
            U[0][k] = (alpha(tau*k)+alpha2(tau*k)*U[1][k]/h) / (alpha1(tau*k)+alpha2(tau*k)/h)
            U[N][k] = (beta(tau*k)+beta2(tau*k)*U[N-1][k]/h) / (beta1(tau*k)+beta2(tau*k)/h)
    elif sigma == 1/2 or sigma == 1:
        for k in range(1,M+1):
            A,B,C,G = [],[],[],[]
            A.append(0)
            B.append(-alpha1(tau*k)-alpha2(tau*k)/h)
            C.append(-alpha2(tau*k)/h)
            G.append(alpha(tau*k))
            tk = tau*k-tau/2 if sigma == 1/2 else tau*k
            for i in range(1,N):
                A.append(sigma*a(h*i,tau*k)/(h**2) - sigma*b(h*i,tau*k)/(2*h))
                B.append(2*sigma*a(h*i,tau*k)/(h**2) - sigma*c(h*i,tau*k)+ 1/tau)
                C.append(sigma*a(h*i,tau*k)/(h**2) + sigma*b(h*i,tau*k)/(2*h))
                G.append(-U[i][k-1]/tau - (1-sigma)*Lh(i,k-1,h,tau,U) - f(h*i,tk))
            A.append(-beta2(tau*k)/h)
            B.append(-beta1(tau*k)-beta2(tau*k)/h)
            C.append(0)
            G.append(beta(tau*k))
            for i in range(N+1): U[i][k] = sweep(N,A,B,C,G)[i]
    return U

def norm(u1,u2):
    mx = 0
    k1,k2,k3,k4 = (len(u2)-1)//5,(len(u2[0])-1)//5,(len(u1)-1)//5,(len(u1[0])-1)//5
    for i in range(6):
        for j in range(6):
            if abs(u1[i*k3][j*k4] - u2[i*k1][j*k2]) > mx:
                mx = abs(u1[i*k3][j*k4] - u2[i*k1][j*k2])
    return mx

def table1(N,M,U):
    k1 = (len(U)-1)//5
    k2 = (len(U[0])-1)//5
    W = grid(5,5)
    Us = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            Us[i][j] = U[i*k1][j*k2]
    v = np.concatenate((np.array([[W[i][0][0] for i in range(len(W))]]).T,rnd(Us,5)),axis=1).T
    trace = go.Table(
        header=dict(values=['x/t']+[str(W[0][i][1]) for i in range(len(W[0]))],
                fill = dict(color=['#e1bccf','#fbcfe5']),
                font = dict(size = 24),
                height = 60),
        cells=dict(values= v.tolist(),
                font = dict(size = 24),
                height = 60,
                fill = dict(color = ['#fbcfe5', 'white']))) 
    data = [trace] 
    plotly.offline.plot(data)    

def table2(U1_lst,U2_lst):
    norm1,norm2 = [norm(U1_lst[i],ex_sol) for i in range(3)],[norm(U1_lst[i],U2_lst[i]) for i in range(3)]
    trace = go.Table(
        header=dict(values=['h','tau','||ex - u||','||u - u2||'],
                fill = dict(color=['#e1bccf','#fbcfe5']),
                font = dict(size = 24),
                height = 60),
        cells=dict(values=[[0.2,0.1,0.05],[0.02,0.005,0.00125],list(map(lambda x: round(x,8),norm1)),list(map(lambda x: round(x,8),norm2))],
                font = dict(size = 24),
                height = 60,
                fill = dict(color = ['#fbcfe5', 'white'])))    
    data = [trace]
    plotly.offline.plot(data)    

def plot_table(N_lst,M,scheme,table_num): #table_num: (0,1,2,3) = (table2,table1(5),table1(10),table1(20))
    if scheme == 2: sigma = 1/2 #change sigma here
    U1_lst = []
    U2_lst = []
    for i in range(len(N_lst)):
        N = N_lst[i]       
        if scheme == 1 or sigma == 0: M = get_M(N,T)
        h,tau = 1/N,T/M
        U = explicit_scheme(N,M) if scheme == 1 else weights_scheme(N,M,sigma)
        U1_lst.append(U)
        if i > 0: U2_lst.append(U)
        if i+1 == table_num: table1(N,M,U)
    N = N_lst[-1]*2       
    if scheme == 1 or sigma == 0: M = get_M(N,T)
    h,tau = 1/N,T/M
    U = explicit_scheme(N,M) if scheme == 1 else weights_scheme(N,M,sigma)
    U2_lst.append(U)
    if table_num == 0: table2(U1_lst,U2_lst)

#explicit scheme:
#plot_table([5,10,20],0,1,1) #N = 5
#plot_table([5,10,20],0,1,2) #N = 10
plot_table([5,10,20],0,1,3) #N = 20
#plot_table([5,10,20],0,1,0) #table for all

#weights scheme:
#plot_table([5,10,20],10,2,1) #N = 5, M = 10
#plot_table([5,10,20],10,2,2) #N = 10, M = 10
#plot_table([5,10,20],10,2,3) #N = 20, M = 10
#plot_table([5,10,20],10,2,0) #table for all
