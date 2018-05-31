import numpy as np
import math

def vector(x): return (np.array([x])).T

def get_eig(A,e):
    n = len(A)
    X = np.identity(n)
    mx = 1000000
    while mx > e: 
        mx,ik,jk = 0,0,0
        for i in range(n):
            for j in range(i + 1,n):
                if abs(A[i][j]) > mx:
                    mx,ik,jk = abs(A[i][j]),i,j
        if mx > e:
            d = math.sqrt((A[ik][ik] - A[jk][jk]) ** 2 + 4 * A[ik][jk] ** 2)
            c = math.sqrt((1 + abs(A[ik][ik] - A[jk][jk]) / d) / 2)
            s = np.sign(A[ik][jk] * (A[ik][ik] - A[jk][jk])) * math.sqrt((1 - abs(A[ik][ik] - A[jk][jk]) / d) / 2)
            V = np.identity(n)
            V[ik][ik],V[jk][jk] = c,c
            V[ik][jk],V[jk][ik] = -s,s
            A2 = np.copy(A)
            for i in range(n):
                if i != ik and i != jk:
                    A[i][ik] = A[ik][i] = c * A2[i][ik] + s * A2[i][jk]
                    A[i][jk] = A[jk][i] = c * A2[i][jk] - s * A2[i][ik]              
            A[ik][ik] = c * c * A2[ik][ik] + 2 * c * s * A2[ik][jk] + s * s * A2[jk][jk]
            A[jk][jk] = s * s * A2[ik][ik] - 2 * c * s * A2[ik][jk] + c * c * A2[jk][jk]
            A[ik][jk] = 0
            A[jk][ik] = 0
            X = np.matmul(X,V)
    X = X.T
    for i in range(len(X)):
        X[i] = X[i] / np.linalg.norm(X[i])
    return np.diagonal(A),X.T

def eig_max_pow(A,e): 
    n = len(A)
    y = vector([1 for i in range(n)])
    eig_mx = [0 for i in range(n)]
    cnt = 0
    while 1:
        cnt+=1
        tmp_y,tmp_mx = y,eig_mx
        y = np.matmul(A,y)
        eig_mx = [y[i]/tmp_y[i] for i in range(n)]
        if abs(eig_mx[0] - tmp_mx[0]) < e and abs(eig_mx[1] - tmp_mx[1]) < e and abs(eig_mx[2] - tmp_mx[2]) < e: break
    return max(eig_mx,key=abs)[0], y/np.linalg.norm(y), cnt

def eig_max_scal(A,e): 
    n = len(A)
    y = vector([1 for i in range(n)])
    eig_mx = 0
    cnt = 0
    while 1:
        cnt+=1
        tmp_y,tmp_mx = y,eig_mx
        y = np.matmul(A,y)
        eig_mx = np.inner(y.T,tmp_y.T)/np.inner(tmp_y.T,tmp_y.T)
        if abs(eig_mx - tmp_mx) < e: break
    return eig_mx[0][0], y / np.linalg.norm(y), cnt
    
def spectrum_bound(A,e):
    eig_A = eig_max_pow(A,e)[0]
    n = len(A)
    B = np.subtract(A,eig_A*np.identity(n))
    eig_B = eig_max_pow(B,e)[0]
    return eig_A + eig_B,eig_max_pow(np.add(A,B),e)[1]

def vilandt(A,eig_0,e):
    n = len(A)
    y = vector([1 for i in range(n)])
    eig_mx = eig_0
    cnt = 0
    while 1:
        cnt+=1
        W = np.subtract(A,eig_mx*np.identity(len(A)))
        y = eig_max_scal(np.linalg.inv(W),e)[1]
        tmp_mx = eig_mx
        eig_mx = 1/eig_max_scal(np.linalg.inv(W),e)[0] + eig_mx
        if abs(eig_mx - tmp_mx) < e: break
    return eig_mx, y, cnt
  
A = np.array([[-0.81417,-0.01937,0.41372],[-0.01937,0.54414,0.00590],[0.41372,0.00590,-0.81445]])
A2 = np.copy(A)
print("A =\r\n",A)
val,vec = get_eig(A2,10e-6)
print ("1) eigenvalues:",val)
print("eigenvectors:\r\n",vec)
#print("check X_inv * A * X:\r\n",np.matmul(np.matmul(np.linalg.inv(vec),A),vec))

eig_p,y_p,cnt_p = eig_max_pow(A,10e-3)
print("\r\n2) max by abs eigenvalue =",eig_p,"\r\neigenvector =\r\n", y_p)
print("neviazka:",np.linalg.norm(np.subtract(np.matmul(A,y_p),eig_p*y_p)))
print("pow iter num =",cnt_p)

eig_s,y_s,cnt_s = eig_max_scal(A,10e-6)
print("\r\n3) max by abs eigenvalue =",eig_s,"\r\neigenvector =\r\n", y_s)
print("neviazka:",np.linalg.norm(np.subtract(np.matmul(A,y_s),eig_s*y_s)))
print("scal iter num =",cnt_s)

eig_op,y_op = spectrum_bound(A,10e-3)
print("\r\n4) spectrum bound =",eig_op,"\r\neigenvector =\r\n",y_op)

eig_v,y_v,cnt_v = vilandt(A,eig_s,10e-3)
print("\r\n5) eigenvalue by vilandt:",eig_v,"\r\neigenvector =\r\n",y_v)
print("neviazka:",np.linalg.norm(np.subtract(np.matmul(A,y_v),eig_v*y_v)))

    
