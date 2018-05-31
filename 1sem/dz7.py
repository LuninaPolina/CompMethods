import math
import numpy.linalg
from scipy.integrate import quad
import numpy
alpha=-0.25 
b=0.5
f=math.cos

def nuka(alpha,k,b): 
    return (b**(alpha+1+ k))/(1+alpha+k)  

def roots2(nu0,nu1,nu2,nu3):
    a1=(nu0*nu3-nu2*nu1)/(nu1**2-nu2*nu0) 
    a2=(nu2**2-nu3*nu1)/(nu1**2-nu2*nu0) 
    d=math.sqrt(a1**2-4*a2)
    x1=(-a1+d)/2 
    x2=(-a1-d)/2
    
    return (x1,x2)

def roots4(nu0,nu1,nu2,nu3,nu4,nu5,nu6,nu7): 
    systemDet=numpy.linalg.det(numpy.array([[nu3,nu2,nu1,nu0],[nu4,nu3,nu2,nu1],[nu5,nu4,nu3,nu2],[nu6,nu5,nu4,nu3]]))
    a1=numpy.linalg.det(numpy.array([[-nu4,nu2,nu1,nu0],[-nu5,nu3,nu2,nu1],[-nu6,nu4,nu3,nu2],[-nu7,nu5,nu4,nu3]]))/systemDet
    a2=numpy.linalg.det(numpy.array([[nu3,-nu4,nu1,nu0],[nu4,-nu5,nu2,nu1],[nu5,-nu6,nu3,nu2],[nu6,-nu7,nu4,nu3]]))/systemDet
    a3=numpy.linalg.det(numpy.array([[nu3,nu2,-nu4,nu0],[nu4,nu3,-nu5,nu1],[nu5,nu4,-nu6,nu2],[nu6,nu5,-nu7,nu3]]))/systemDet
    a4=numpy.linalg.det(numpy.array([[nu3,nu2,nu1,-nu4],[nu4,nu3,nu2,-nu5],[nu5,nu4,nu3,-nu6],[nu6,nu5,nu4,-nu7]]))/systemDet
    (x1,x2,x3,x4) = numpy.roots([1,a1,a2,a3,a4])
    
    return (x1,x2,x3,x4)

x1,x2=roots2(nuka(alpha,0,b),nuka(alpha,1,b),nuka(alpha,2,b),nuka(alpha,3,b)) 
print(x1,x2)
y1,y2,y3,y4=roots4(nuka(alpha,0,b),nuka(alpha,1,b),nuka(alpha,2,b),nuka(alpha,3,b),nuka(alpha,4,b),nuka(alpha,5,b),nuka(alpha,6,b),nuka(alpha,7,b))
def A1(x):
    return x**alpha*(x-x2)/(x1-x2)
def A2(x):
    return x**alpha*(x-x1)/(x2-x1)
def B1(x):
    return x**alpha*(x-y2)*(x-y3)*(x-y4)/(y1-y2)/(y1-y3)/(y1-y4)
def B2(x):
    return x**alpha*(x-y1)*(x-y3)*(x-y4)/(y2-y1)/(y2-y3)/(y2-y4)
def B3(x):
    return x**alpha*(x-y2)*(x-y1)*(x-y4)/(y3-y2)/(y3-y1)/(y3-y4)
def B4(x):
    return x**alpha*(x-y2)*(x-y3)*(x-y1)/(y4-y2)/(y4-y3)/(y4-y1)
i1=f(x1)*(quad(A1,0,b))[0] + f(x2)*(quad(A2,0,b))[0]
i2=f(y1)*(quad(B1,0,b))[0] + f(y2)*(quad(B2,0,b))[0]+ f(y3)*(quad(B3,0,b))[0] + f(y4)*(quad(B4,0,b))[0]
print("n = 2; integral = ",i1)
print (y1,y2,y3,y4)
print("n = 4; integral = ",i2)






 

