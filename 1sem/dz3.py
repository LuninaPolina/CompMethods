xArr = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
yArr = [1.623250, 1.664792, 1.701977, 1.734832, 1.763404, 1.787764, 1.808002, 1.824230, 1.836580, 1.845201]
x1 = 0.159884
x2 = 0.801404
x3 = 0.207736
y = 1.792942
x0 = xArr[1]
y0 = yArr[1]
xn = xArr[9]
yn = yArr[9]

def createTable (xArr, yArr):
    table = []
    table.append(xArr)
    table.append(yArr)
    tmp = yArr
    cnt = len(yArr)-1
    while not (all(el == 0 for el in table[len(table) - 1])):
        table.append([tmp[i+1]-tmp[i] for i in range(cnt-1)] + [0 for i in range(len(yArr)-1-cnt)])
        tmp = table[len(table)-1]
        cnt-=1
    zip(table)
    return table
table = createTable(xArr,yArr)
t1 = (x1-xArr[1])/0.1
p1 = y0 + t1*table[1][2] + t1*(t1-1)*table[1][3]/2 + t1*(t1-1)*(t1-2)*table[1][4]/6 + t1*(t1-1)*(t1-2)*(t1-3)*table[1][5]/24
t2 = (x2-xn)/0.1
p2 = yn + t2*table[9][2] + t2*(t2-1)*table[8][3]/2 + t2*(t2-1)*(t2-2)*table[7][4]/6 + t2*(t2-1)*(t2-2)*(t2-3)*table[6][5]/24
print(p1,p2)


    
