import random
import numpy as np
from numpy import sin, cos
from itertools import product

import matplotlib.pyplot as plt

#IMPORTANT
#np.seterr(all='raise')

# Dot product of unit vector spins s_1 * s_2
def dot_Product(angle_set1, angle_set2):

    theta1, phi1 = angle_set1[0], angle_set1[1]
    theta2, phi2 = angle_set2[0], angle_set2[1]
    dot_product = sin(theta1)*sin(theta2)*cos(phi1)*cos(phi2)+sin(theta1)*sin(theta2)*sin(phi1)*sin(phi2)+cos(theta1)*cos(theta2)

    return dot_product

def max_value_of_Matrix(matrix):

    m = []
    for i in matrix:
        for j in i:
            m.append(j)
    return max(m)

def IsingMatrix(J):

    s_vals = [1,-1]
    hamiltonian = []

    for i in range(2):
        for j in range(2):
            hamiltonian.append(J*(s_vals[i]*s_vals[j]))

    hamiltonian_max = np.amax(hamiltonian)

    Transfer = []
    num = 0
    for i in range(2):
        row = []
        for j in range(2):
            row.append(np.exp(hamiltonian[num]-hamiltonian_max))
            num = num+1
        Transfer.append(row)

    return np.array(Transfer,dtype=np.longdouble)
    
            
def TransferMatrixCreator(J):

    theta_list = np.radians(np.linspace(0,180,10)) #Theta(in degrees) (0,180,10)
    phi_list = np.radians(np.linspace(0,342,20)) #Phi(in degrees) (0,342,20)
    angle_set = list(product(theta_list,phi_list)) #Kartezyen Çarpımı

    hamiltonian = []

    for i in range(len(angle_set)):
        for j in range(len(angle_set)):
            hamiltonian.append(J*dot_Product(angle_set[i],angle_set[j]))

    hamiltonian_max = np.amax(hamiltonian)

    Transfer = []
    num = 0
    for i in range(len(angle_set)):
        row = []
        for j in range(len(angle_set)):
            row.append(np.exp(hamiltonian[num]-hamiltonian_max))
            num = num+1
        Transfer.append(row)

    return np.array(Transfer,dtype=np.longdouble)

def TransferMatrixDec(matrix1, matrix2):
    
    Transfer = []
    Transfer = np.matmul(matrix1, matrix2)
    Transfer = Transfer * (1/max_value_of_Matrix(Transfer))
    return Transfer

def RG_Flow(RG_step, J_initial):

    Flow_TM = []
    tm = TransferMatrixCreator(J_initial)
    #tm = IsingMatrix(J_initial)
    
    print('Mean value of Matrix: {}'.format(np.mean(tm)))
    Flow_TM.append(tm)
    
    for i in range(RG_step):

        tm_transformed = TransferMatrixDec(tm,tm)
        tm = tm_transformed
        Flow_TM.append(tm)
        print('RG STEP NO : {}'.format(i+1))
        print('Mean value of Matrix: {}'.format(np.mean(tm)))
        print('Min value of Matrix: {}'.format(np.amin(tm.flatten())))
        if np.mean(tm) == 1:
            break
        else:
            pass
    
    return Flow_TM
