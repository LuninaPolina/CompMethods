import pylab
import math
def f(x,y):
    return -(y+x)/(y+x+1)
y0=0
a,b=0,1
tmp=[]
def runge(h,xArr):
    yArr=[y0]
    for i in range(1,5):
        k1=h*f(xArr[i-1],yArr[i-1])
        k2=h*f(xArr[i-1]+h/2,yArr[i-1]+k1/2)
        k3=h*f(xArr[i-1]+h/2,yArr[i-1]+k2/2)
        k4=h*f(xArr[i-1]+h,yArr[i-1]+k3)
        yArr.append(yArr[i-1]+(k1+2*k2+2*k3+k4)/6)
    return yArr
def createTable (xArr, yArr,h,n):
    table = []
    table.append(xArr)
    table.append(yArr)
    table.append([])
    for i in range(len(yArr)): table[2].append(h*f(xArr[i],yArr[i]))
    k=4
    for i in range(2,6):
        table.append([table[i][j+1]-table[i][j] for j in range(k)])
        k-=1
    for i in range(5,n+1):
        table[1].append(table[1][i-1]+table[2][i-1]+table[3][i-2]/2+5*table[4][i-3]/12+3*table[5][i-4]/8+251*table[6][i-5]/720)
        table[2].append(h*f(xArr[i],yArr[i]))
        for j in range(3,7):
            if i+3-j<len(table[j-1]): table[j].append(table[j-1][i+3-j]-table[j-1][i+2-j])     
    return table
def eu(n):
    h=(b-a)/n
    xArr=[]
    for i in range(n+1):
        xArr.append(a+i*h)
    yArr=[y0]
    for i in range(1,n+1):
        yArr.append(yArr[i-1]+h*f(xArr[i-1],yArr[i-1]))
    print("n = ",n)
    for i in range(len(xArr)):
        print(i,"  ",round(xArr[i],3),"  ",round(yArr[i],3))
    print("*********************")
    if n==10: tmp.append(yArr[10])
    pylab.plot(xArr,yArr)
def ad(n):
    h=(b-a)/n
    xArr=[]
    for i in range(n+1):
        xArr.append(a+i*h)
    yArr=runge(h,xArr)
    table=createTable(xArr,yArr,h,n)
    print("n = ",n)
    for i in range(len(xArr)):
        print(i,"  ",round(xArr[i],3),"  ",round(yArr[i],3))
    print("*********************")
    if n==10: tmp.append(table[1][10])
    pylab.plot(xArr,table[1])    

eu(10)
eu(20)
eu(5)
pylab.title("Euler")
pylab.show()
ad(10)
ad(20)
ad(5)
pylab.title("Adams")
pylab.show()

A=0.1*math.exp(1)
print("|ya-ye| = ", abs(tmp[0]-tmp[1]))
print("A = ",A)
