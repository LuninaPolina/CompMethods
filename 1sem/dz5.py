import math
def f(x):
    return 1/(1.28-math.exp(1)**(-x-0.01))
a,b=0,0.4
mx=[1.07167/(0.226524**3),1.07167/(0.226524**3),10.002/(0.226524**5)]
def priamougolnik (n):
    s=0
    h=(b-a)/n
    for i in range (n):
        s+=f(a+(h/2)*(1+2*i))
    return s*h
def trapetsia (n):
    s=0
    h=(b-a)/n
    for i in range (1,n):
        s+=f(a+i*h)
    return h*(s+(f(a)+f(b))/2)
def simpson (n):
    h=(b-a)/n
    s=4*(f(a+h/2))    
    for i in range (1,n):
        s=s+2*(f(a+i*h)) + 4*(f(a+h/2*(2*i+1)))
    return h/6*(s+f(a)+f(b))
def runge(j8,j16,k):
    return (j8 - j16*2**k)/(1-2**k)
def r(i):
    h=(b-a)/8
    if i==0: t=(h**2)*(b-a)/24
    elif i==1: t=-(h**2)*(b-a)/12
    else: t=-(h**4)*(b-a)/2880
    return abs(mx[i]*t)
j8=priamougolnik(8)
j16=priamougolnik(16)
print("rectangle:", "\r\nJ8 = ",j8,"\r\nJ16 = ",j16,"\r\nrunge = ",runge(j8,j16,2),"\r\nRn = ",r(0))
j8=trapetsia(8)
j16=trapetsia(16)
print("trapeze:", "\r\nJ8 = ",j8,"\r\nJ16 = ",j16,"\r\nrunge = ",runge(j8,j16,2),"\r\nRn = ",r(1))
j8=simpson(8)
j16=simpson(16)
print("simpson:", "\r\nJ8 = ",j8,"\r\nJ16 = ",j16,"\r\nrunge = ",runge(j8,j16,4),"\r\nRn = ",r(2))
print("********************")
j16=runge(j8,j16,4)
j8=priamougolnik(8)
print("|rungesimpson-rectangle|:", abs(j8-j16))
j8=trapetsia(8)
print("|rungesimpson-trapeze|:", abs(j8-j16))
j8=simpson(8)
print("|rungesimpson-simpson|:", abs(j8-j16))


