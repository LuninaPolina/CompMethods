import math
import numpy as np
import scipy.integrate as integrate
import scipy.special as special


def f(x): return x + 0.1
def H(x,y): return -0.1*(math.sin(x*(0.5+y**2)))
a,b = 0,1

def alpha(i,x):
    return -0.1*((-1)**i)*(x**(2*i+1))/math.factorial(2*i+1)
def beta(i,y):
    return (0.5+y**2)**(2*i+1)


def Hn(x,y,n):
    summa = 0
    for i in range(n):
        summa = summa + alpha(i,x)*beta(i,y)
    return summa

def kron(i,j):
    if i == j: return 1
    else: return 0

def G(x,y,n,D):
    summa = 0
    for i in range(n):
        for j in range(n):
            summa = summa + D[i][j]*alpha(i,x)*beta(j,y)
    return summa           

def un(n):
    A = np.zeros((n, n))
    B = np.zeros((n,1))
    for i in range(n):
        for j in range(n):
            A[i][j] = kron(i,j) - integrate.quad(lambda y: beta(i,y)*alpha(j,y),a,b)[0]
            B[i][0] = integrate.quad(lambda y: beta(i,y)*f(y),a,b,)[0]
    C = np.linalg.solve(A,B)
    D = np.linalg.inv(A)
    B,nu = 0,0
    for x in range(101):
        tmp = integrate.quad(lambda y: abs(G(a+(x/100)*b,y,n,D)),a,b)[0] 
        if tmp > B: B = tmp
        tmp = integrate.quad(lambda y: abs(H(a+(x/100)*b,y)-Hn(a+(x/100)*b,y,n)),a,b)[0]
        if tmp > nu: nu = tmp
    return lambda x: f(x) + integrate.quad(lambda y: G(x,y,n,D)*f(y),a,b)[0],B,nu

def estimate(u,B,nu):
    norm_u = 0
    for x in range(101):
        if abs(u(a+(x/100)*b)) > norm_u: norm_u = abs(u(a+(x/100)*b))
    return (1+B)*nu*norm_u/(1-(1+B)*nu)

def calculate(x,A,X,z):
    summa = 0
    for i in range(len(A)):
        summa = summa + A[i]*H(x,X[i])*z[i]
    return summa[0] + f(x)

def solve(n): #n - number of pieces
    h=(b-a)/n
    X = [a+i*h for i in range(n+1)]
    A = [h for i in range(n+1)]
    A[0],A[n] = h/2,h/2
    D = np.zeros((n+1, n+1))
    g = np.zeros((n+1, 1))
    for i in range(n+1):
        g[i] = f(X[i])
        for j in range(n+1):
            D[i][j] = kron(i,j) - A[j]*H(X[i],X[j])
    z = np.linalg.solve(D,g)
    return lambda x: calculate(x,A,X,z)
    
def task_1():
    u3,B3,nu3 = un(3)
    print("x           a           (a+b)/2        b")
    print("u3(x)      ",round(u3(a),8),"       ",round(u3((a+b)/2),8),"   ",round(u3(b),8))
    u4,B4,nu4 = un(4)
    print("u4(x)      ",round(u4(a),8),"       ",round(u4((a+b)/2),8),"   ",round(u4(b),8))
    print("delta      ",max(abs(u4(a)-u3(a)),abs(u4((a+b)/2)-u3((a+b)/2)),abs(u4(b)-u3(b))))
    print("aposterior ",estimate(u3,B3,nu3))

def task_2():
    print("x           a           (a+b)/2        b")
    u2 = solve(2)
    print("u2(x)      ",round(u2(a),8),"       ",round(u2((a+b)/2),8),"   ",round(u2(b),8))
    u4 = solve(4)
    print("u4(x)      ",round(u4(a),8),"       ",round(u4((a+b)/2),8),"     ",round(u4(b),8))
       
    n = 8
    while max(abs(u4(a)-u2(a)),abs(u4((a+b)/2)-u2((a+b)/2)),abs(u4(b)-u2(b))) > 10e-5:
       u2 = u4
       u4 = solve(n)
       print("u"+str(n)+"(x)     ",round(u4(a),8),"       ",round(u4((a+b)/2),8),"   ",round(u4(b),8))
       n *=2
    print("u"+str(n)+"(x)-u"+str(n-1)+"(x)   ",max(abs(u4(a)-u2(a)),abs(u4((a+b)/2)-u2((a+b)/2)),abs(u4(b)-u2(b))))

num = input("Task number: ")
if num == "1": task_1()
if num == "2": task_2()    

            


