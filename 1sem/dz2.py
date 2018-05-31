import numpy
import math

def f(p,k):
    (x,y) = (p[0],p[1])
    return math.tan(x-y+k)-x*y

def g(p,a):
    (x,y) = (p[0],p[1])
    return a*x*x+2*y*y-1

def cosinus (x,y,k):
    return (math.cos(x-y+k))*(math.cos(x-y+k)) 

def delta(p,k,a):
    (x,y) = (p[0],p[1])
    return 4*y/cosinus(x,y,k)-4*y*y+2*a*x/cosinus(x,y,k)+2*a*x*x

def deltax(p,k,a):
    (x,y) = (p[0],p[1])
    return -f(p,k)*4*y-g(p,a)/cosinus(x,y,k)-g(p,a)*x

def deltay(p,k,a):
    (x,y) = (p[0],p[1])
    return -g(p,a)/cosinus(x,y,k)+g(p,a)*y+f(p,k)*2*a*x

def d (lx,rx,ly,ry,n):
    points = []
    for i in range (n+1):
        points.append([])
        for j in range (n+1):
            points[i].append(((lx+(rx-lx)/n*i,ly+(ry-ly)/n*j)))
    return points

def find (k,a,lx,rx,ly,ry):
    points = d(lx,rx,ly,ry,1000)
    tmp = []
    for i in range(0, len(points) - 2):
        for j in range(0, len(points[0]) - 2):
            if (f(points[i][j],k)*f(points[i+1][j+1],k) < 0 and g(points[i][j],a)*g(points[i+1][j+1],a) < 0) or (f(points[i][j],k)*f(points[i+1][j],k) < 0 and g(points[i][j],a)*g(points[i+1][j],a) < 0) or (f(points[i][j],k)*f(points[i][j+1],k) < 0 and g(points[i][j],a)*g(points[i][j+1],a) < 0):
                tmp.append((points[i][j],points[i+1][j+1]))
    for el in tmp:
        for el2 in tmp:
            if abs(el2[0][0]-el[0][0]) < 0.1 and abs(el2[0][1]-el[0][1]) < 0.1 and el != el2:
                tmp.remove(el2)
    for el in tmp:
        cnt = 1
        if abs(f(el[0],k)) <= 0.1 and abs(g(el[0],a)) <= 0.1:
            p = el[0]
            print(cnt,p,f(p, k),g(p, a))
            while (abs(f(p,k)) > 0.00001 and abs(g(p,a)) > 0.00001) or (abs(deltax(p,k,a)/delta(p,k,a)) > 0.0000000001 and abs(deltay(p,k,a)/delta(p,k,a)) > 0.0000000001):               
                cnt +=1
                p = (p[0]+deltax(p,k,a)/delta(p,k,a),p[1]+deltay(p,k,a)/delta(p,k,a))
                print(cnt,p,f(p, k),g(p, a))
            print ("-----------\r\n")
        else:
            find(k,a,el[0][0],el[1][0],el[0][1],el[1][1])
achki = [0.5,0.6,0.7,0.8,0.9,1.0]
kachki = [0.2,0.1,0.0,-0.1,-0.2]
for i in range (0,6):
    for j in range (0,5):
        print ("a = ",achki[i], "k = ",kachki[j])
        find (kachki[j],achki[i],-2,2,-2,2)
