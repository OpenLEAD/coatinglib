from numpy import *
Y = load('trajectory/Y.npz')
Y = Y['array']

y=[]
tempy=[]
comp = Y[0][0]
for i in Y[0]:
    if abs(sqrt(i[0]**2+i[1]**2+i[2]**2)-sqrt(comp[0]**2+comp[1]**2+comp[2]**2))<=0.00001:
        tempy.append(i)
    else:
        comp = i
        y.append(tempy)
        tempy=[]

savez_compressed('trajectory/'+'Y2.npz', array=y)         
