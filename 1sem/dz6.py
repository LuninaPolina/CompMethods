import math

x1=7/30 + (2*math.sqrt(14/11))/15
x2=7/30 - (2*math.sqrt(14/11))/15
J = math.cos(x1)*(0.169887-0.792805*x2)/(x1-x2) + math.cos(x2)*(-0.169887+0.792805*x1)/(x1-x2)
print(J)
