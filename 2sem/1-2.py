import numpy as np
import math

def vector(x): return (np.array([x])).T

def diagonal(A):
    D = np.copy(A)
    for i in range(len(A)):
        for j in range(len(A[0])):
            if i != j: D[i][j] = 0
    return D
    
def inverse_by_cofactors(A): #only 2x2 matrix
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    inv = np.zeros([2,2])
    for i in range(2):
        for j in range(2):
            inv[i][j] = int(((-1) ** (i + j)) * A[1 - i][1 - j] / det)
    return inv

def norm(A): #vector(x) if x is vector
    rows = []
    for i in range(len(A)): rows.append(sum(np.absolute(A[i])))
    return max(rows)

def forward(A,b):
    n = len(A)
    m = len(b.T)
    A_b = np.concatenate((A,b),axis=1)
    for k in range(n):
        if k != n - 1: choose_main_el(A_b,k) 
        tmp = A_b[k][k]
        if abs(tmp) < 0.0001:
            print('Too small main element!')
        for j in range(k + 1,n + m): #from k if all
            A_b[k][j] = A_b[k][j] / tmp
        for i in range(k + 1,n):
            tmp = A_b[i][k]
            for j in range(k + 1,n + m): #from k if all
                A_b[i][j] = A_b[i][j] - A_b[k][j] * tmp
    return A_b

def backward(A_b):
    n = len(A_b)
    x = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    for l in range(len(A_b[0])-len(A_b)):
        x[l][n - 1] = A_b[n - 1][n]
        for i in range(n - 2,-1,-1):
            s = 0
            for j in range(i + 1,n):
                s = s + A_b[i][j] * x[l][j]
            x[l][i] = A_b[i][n] - s
    return x

def choose_main_el(A_b,k):
    mx = 0
    for i in range(k,len(A_b)):
        if A_b[i][k] >= mx:
            mx = A_b[i][k]
            ind = i
    tmp = np.copy(A_b[ind])
    A_b[ind] = A_b[k]
    A_b[k] = tmp

def gauss(A,b):
    return backward(forward(np.copy(A),np.copy(b)))

def inverse_by_gauss(A):
    n = len(A)
    res = np.zeros([n,n])
    b = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    res = gauss(A,b)
    return np.transpose(res)

def aprior(H,g):
    res = norm(H) * norm(g) / (1 - norm(H))
    k = 1
    while res >= 0.001:
        k += 1
        res = (norm(H) ** k) * norm(g) / (1 - norm(H))
    return k, res

def gershgorin(A):
    left,right = [],[]
    for i in range(len(A)):
        a = A[i][i]
        s = 0
        for j in range(len(A[0])):
            if i != j: s += abs(A[i][j])
        left.append(a - s)
        right.append(a + s)
    return min(left),max(right)

def simple_iter(H,g,k):
    x = vector([0 for i in range(len(H))])
    for i in range(1,k + 1):
           x = np.add(np.matmul(H,x),g)
    return x

def lusternik(H,g,xk,k):
    p = max(np.absolute(np.linalg.eig(H)[0]))
    x = simple_iter(H,g,k - 1) + (xk - simple_iter(H,g,k - 1))/(1 - p)
    return x

def seidel(H,g,k):
    n = len(H)
    HL = np.array([[0,0,0],[H[1][0],0,0],[H[2][0],H[2][1],0]]) #fix
    HR = np.array([[H[0][0],H[0][1],H[0][2]],[0,H[1][1],H[1][2]],[0,0,H[2][2]]])
    x = vector([0 for i in range(n)])
    for k in range(1,k + 1):
        tmp = np.linalg.inv(np.subtract(np.identity(n),HL))
        x = np.add(np.matmul(np.matmul(tmp,HR),x),np.matmul(tmp,g))
    return x   

def seidel_by_formula(H,g,k):
    n = len(H)
    x = vector([0.0 for i in range(n)])
    for l in range(1,k + 1):
        for i in range(n):
            s1,s2 = 0,0
            for j in range(i):
                s1 = s1 + H[i][j] * x[j]
            for j in range(i,n):
                s2 = s2 + H[i][j]* x[j]
            x[i][0] = s1 + s2 + g[i][0]
    return x

def cheb_iter(A,b,p):
    n = len(A)
    m,M = gershgorin(A)
    l = math.ceil(math.log(p,2))
    x = vector([0 for i in range(n)])
    tetas = []
    tetas.append([])
    tetas[0].append(1)
    for i in range(1,l + 1): 
        tetas.append([0 for k in range(2 ** i)])
        for j in range(0,2 ** (i - 1) ):
            tetas[i][2 * j] = tetas[i - 1][j]
            tetas[i][2 * j + 1] = 4 * 2 ** (i - 1) - tetas[i][2 * j]
    stop = p
    for k in range(stop): 
        t = math.cos(math.pi * (tetas[-1][k]) / 2 / p)
        tau = 2 / (M + m - (M - m) * t)
        x = x + tau * (b - np.matmul(A,x))
    return x

def task_1_1():
    A = np.array([[1,0.99],[0.99,0.98]])
    b = vector([1.99,1.97])
    print("A_b =\r\n",np.concatenate((A,b),axis=1))
    db = vector([0.01,0.03])
    print("b2 =",np.add(db,b).T[0])
    x = np.linalg.solve(A,b)
    dx = np.subtract(np.linalg.solve(A,np.add(b,db)),x)
    mu = norm(A) * norm(inverse_by_cofactors(A))
    delta_x = norm(dx) / norm(x)
    delta_b = norm(db) / norm(b)
    print("mu(A) =",mu)
    print("deltaB =",delta_b)
    print("deltaX =",delta_x)
    print(delta_x,"<=",mu * delta_b,"=> deltaX <= mu(A) * deltaB")

def task_1_2(): #var 2
    A = np.array([[7.35272,0.88255,-2.270052],[0.88255,5.58351,0.528167],[-2.27005,0.528167,4.430329]])
    b = vector([1,0,0])
    A_b = np.concatenate((A,b),axis=1)
    print("A_b =\r\n",A_b)
    x = gauss(A,b)
    print("x =",x[0])
    print("R =",np.subtract(b,np.matmul(A,x[0])).T[0])
    A_inv = inverse_by_gauss(A)
    print("A_inv =\r\n",A_inv)

def task_2():
    A = np.array([[12.785723,1.534675,-3.947418],[1.534675,9.709232,0.918435],[-3.947418,0.918435,7.703946]])
    b = vector([9.60565,7.30777,4.21575])
    A_b = np.concatenate((A,b),axis=1)
    x = gauss(A,b)
    print("A_b =\r\n",A_b)
    print("\r\nx =",(x.T)[0])
    D = diagonal(A)
    g = np.matmul(inverse_by_gauss(D),b)
    H = np.subtract(np.identity(len(A)),np.matmul(inverse_by_gauss(D),A))

    print("\r\nH =\r\n",H)
    print("norm(H) =",norm(H))
    
    print("\r\nk =", aprior(H,g)[0], "aprior:",aprior(H,g)[1])
    xk = simple_iter(H,g,1)
    k = 1
    while norm(np.subtract(xk,x)) >= 0.001:
            k += 1
            xk = simple_iter(H,g,k)
    print("\r\nx_iter =",(xk.T)[0],"k =",k)
    print("||x_iter - x|| = ",norm(np.subtract(x,xk)))
    print("aprior:", (norm(H) ** k) * norm(g) / (1 - norm(H)))
    print("aposterior:",  norm(H) * norm(np.subtract(xk,simple_iter(H,g,k - 1))) / (1 - norm(H)))
    xl = lusternik(H,g,xk,k)
    print("lusternik:",(xl.T)[0])
    print("||x_lust - x|| =",norm(np.subtract(x,xl)))

    xs = seidel(H,g,1)
    k2 = 1
    while norm(np.subtract(xs,x)) >= 0.001:
            k2 += 1
            xs = seidel(H,g,k2)
    xs2 = seidel_by_formula(H,g,k2)
    print("\r\nseidel:",xs.T[0], "k =",k2)
    print("seidel by formula:",xs2.T[0], "k =",k2)
    print("||x_seid - x|| =",norm(np.subtract(x,xs)))

    xc = cheb_iter(A,b,2)
    k3 = 2
    while norm(np.subtract(xc,x)) >= 0.001:
            k3 *= 2
            xc = cheb_iter(A,b,k3)
    print("\r\nx_cheb:",xc.T[0],"k =",k3)
    print("||x_cheb - x|| =",norm(np.subtract(x,xc)))
    
num = input("Task number: ")
if num == "1.1": task_1_1()
if num == "1.2": task_1_2()
if num == "2": task_2()



