import numpy as np
import matplotlib.pyplot as plt

def g(J,beta):
    
    return 4*np.pi*np.sinh(J*beta)*(1/(J*beta))

def fEnergy(J,beta):

    return -1*np.log(g(J,beta))

def der(J,beta):

    eps = 10e-8
    diff = fEnergy(J,beta+eps)-fEnergy(J,beta)
    return diff/eps


x = np.linspace(0.1,100,5000)
plt.plot(1/x,-1*(x**2)*der(1,der(1,x)))
plt.show()
