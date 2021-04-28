import random
import numpy as np
from numpy import sin, cos
from itertools import product

from scipy.special import sph_harm, spherical_jn

import pandas as pd
import matplotlib.pyplot as plt

def exp_of_dot_Product_sph(angle_set1, angle_set2, J):

    theta1, phi1 = angle_set1[0], angle_set1[1]
    theta2, phi2 = angle_set2[0], angle_set2[1]

    x = 0
    for l in range(4):
        for m in np.linspace(-l,l,2*l+1):        
            x += 4*np.pi*(1j**l)*spherical_jn(l,-1*1j*J)*np.conj(sph_harm(m,l,phi1,theta1))*sph_harm(m,l,phi2,theta2)

    return x
    
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

    theta_list = np.radians(np.linspace(0,180,180)) #Theta(in degrees) (0,180,10)
    phi_list = np.radians(np.linspace(0,342,20)) #Phi(in degrees) (0,342,20)
    angle_set = list(product(theta_list,phi_list)) #Kartezyen Çarpımı

    hamiltonian = []

    for i in range(len(angle_set)):
        for j in range(len(angle_set)):
            
            hamiltonian.append(J*dot_Product(angle_set[j],angle_set[i]))
            #hamiltonian.append(exp_of_dot_Product_sph(angle_set[i],angle_set[j],J))
            
    hamiltonian_max = np.amax(hamiltonian)

    Transfer = []
    num = 0
    for i in range(len(angle_set)):
        row = []
        for j in range(len(angle_set)):
            row.append(np.exp(hamiltonian[num]-hamiltonian_max))
            #row.append(hamiltonian[num]/hamiltonian_max)
            num = num+1
        Transfer.append(row)

    return np.array(Transfer)

def TransferMatrixDec(matrix1, matrix2):
    
    Transfer = []

    Transfer = np.dot(matrix1, matrix2)

    #if np.mean(Transfer) == np.inf:
    #    Transfer = np.dot(matrix1/np.amax(matrix1),matrix2/np.amax(matrix2))
        
    Transfer = Transfer * (1/max_value_of_Matrix(Transfer))
    Transfer = np.around(Transfer, decimals=8)
    
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


def flow_to_excel(flow):

    df_list = []
    file_name = 'matrice.xlsx'
    writer = pd.ExcelWriter(file_name)
    for i in range(len(flow)):
        df_list.append(pd.DataFrame(np.real(flow[i])))
        df_list[i].to_excel(writer, sheet_name='RG_NO_{}'.format(i), float_format='%.5f')
    writer.save()

    return True

def flow_reduce(Flow_TM):

    reduced_flow = []
    for i in range(len(Flow_TM)):
        reduced_matrix = []
        for j in range(10):
            row = []
            for k in range(10):
                row.append(Flow_TM[i][j*20][k*20])
            reduced_matrix.append(row)
        reduced_flow.append(reduced_matrix)
        reduced_flow[i] = np.array(reduced_flow[i]).reshape(10,10)

    return reduced_flow



    
        
