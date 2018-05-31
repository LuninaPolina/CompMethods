import math
import scipy.interpolate
import numpy
import pylab
from scipy.interpolate import spline

nodes=[math.pi/5,math.pi/3,3*math.pi/4,math.pi]
def createNodes(a,b,n):
    nodes=[]
    for i in range(n):
      nodes.append((b-a)*i/n+a)
    return nodes
def kostil(x,n):
    if n==4: return abs(math.cos(x))
    else: return abs(-2097152*math.cos(8*x))
points=[]
a,b=0,math.pi
for i in range(21):
    points.append((b-a)*i/20+a)
def w(x,n):
    res=1
    for i in range(n):
        res=res*(x-nodes[i])
    return res
def createTable (f,nodes,n,title):
    table=[]
    table.append(points)
    table.append(list(map(lambda x: f(x),points)))
    values_in_nodes=list(map(lambda x: f(x),nodes))
    table.append(list(map(lambda x: numpy.polyval(scipy.interpolate.lagrange(nodes,values_in_nodes),x),points)))
    print(scipy.interpolate.lagrange(nodes,values_in_nodes))
    abses=[]
    for i in range(21):
        abses.append(abs((table[1][i]-table[2][i])))
    table.append(abses)
    mx=80000#max(list(map(lambda x: kostil(x,n),nodes)))
    a=list(map(lambda x: abs(w(x,n))*mx/math.factorial(n),points))
    table.append(a)

    print(title)
    print("xi             f(xi)          Ln(xi,f)       |Ln-f|         A")
    for i in range(21):
        row = ""
        for j in range(len(table)):
            k=15-len(str(round(table[j][i],5)))
            row=row+str(round(table[j][i],5))+ " "*k
        print(row)
    xnew = numpy.linspace(min(table[0]),max(table[0]),300)
    pylab.plot(xnew,spline(table[0],table[1],xnew))
    pylab.plot(xnew,spline(table[0],table[2],xnew))
    pylab.ylim((-1,1))
    pylab.title(title)
    pylab.show()

'''
createTable(math.sin,nodes,4,"sin x, normal")
createTable(math.sin,createNodes(0,math.pi/2,4),4,"sin x, left")
createTable(math.sin,createNodes(math.pi/2,math.pi,4),4, "sin x, right")
createTable(math.sin,createNodes(math.pi/4,3*math.pi/4,4),4,"sin x, middle")
createTable(math.sin,createNodes(0,math.pi,8),4,"sin x, double")    
'''
points=[]
a,b=-math.pi/2,math.pi/2
for i in range(21):
    points.append((b-a)*i/20+a)
def sin8(x):
    return math.sin(10*x)
nodes=[-math.pi/3,-math.pi/8,math.pi/6,math.pi/4,math.pi/2]#[-0.97,-0.69,-0.21,0.33,0.74,0.95]

createTable(sin8,nodes,5,"sin 8x, normal")
createTable(sin8,createNodes(a,b/2,5),5,"sin 8x, left")
createTable(sin8,createNodes(a/2,b,5),5,"sin 8x, right")
createTable(sin8,createNodes(-1/2,1/2,5),5,"sin 8x, middle")
createTable(sin8,createNodes(a,b,10),5,"sin 8x, double")
