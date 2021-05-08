import numpy as np
from scipy.special import spherical_jn
import matplotlib.pyplot as plt

def g(l,J):
    
    return 4*np.pi*(1j**l)*spherical_jn(l,-1j*J)




x = np.linspace(0,10,500)
plt.plot(x,-1*np.log(g(0,x)))
plt.show()

