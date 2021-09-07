from random import uniform
import matplotlib.pyplot as plt
import numpy as np
L=[]
V=[]
Vmoy=[0,0]
pop=1000
xmin, ymin, xmax, ymax = -1.1, -1.1, 1.1, 1.1
plt.xlim(xmin, xmax) # Limites des axes du repère
plt.ylim(ymin, ymax)
for i in range(pop):
    O=uniform(-np.pi/6,np.pi/6)
    L.append(O)
    V.append([np.cos(O),np.sin(O)])
    Vmoy[0]+=np.cos(O)/pop
    Vmoy[1]+=np.sin(O)/pop
    plt.plot(np.cos(O),np.sin(O),'bo')
plt.plot(Vmoy[0],Vmoy[1],'ro')
plt.show()
print(np.sqrt((Vmoy[0])**2+(Vmoy[1])**2))