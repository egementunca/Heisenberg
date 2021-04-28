import numpy as np
from scipy.special import spherical_jn
import matplotlib.pyplot as plt

def lamda(l,J):
    
    return 4*np.pi*(1j**l)*spherical_jn(l,-1j*J)

def eigen(J):
    
    kappa_steps = []
    for l in range(100):
        x = (lamda(l,J)/lamda(0,J))**2
        kappa_steps.append(x)
    return kappa_steps, np.sum(kappa_steps)


