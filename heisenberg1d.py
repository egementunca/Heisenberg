import random
import numpy as np
from numpy import sin, cos
from itertools import product

from scipy.special import sph_harm, spherical_jn, eval_legendre

import pandas as pd
import matplotlib.pyplot as plt

#IMPORTANT
#np.seterr(all='raise')

# Dot product of unit vector spins s_1 * s_2
def exp_of_dot_Product_sph(angle_set1, angle_set2, J):

    theta1, phi1 = angle_set1[0], angle_set1[1]
    theta2, phi2 = angle_set2[0], angle_set2[1]

    x = 0
    for l in range(40):
        for m in np.linspace(-l,l,2*l+1):
            x += 4*np.pi*(1j**l)*spherical_jn(l,-1*1j*J)*np.conj(sph_harm(m,l,phi1,theta1))*sph_harm(m,l,phi2,theta2)
    
    return x

def exp_of_dot_Product_legendre(angle_set1, angle_set2, J):

    theta1, phi1 = angle_set1[0], angle_set1[1]
    theta2, phi2 = angle_set2[0], angle_set2[1]

    gamma = cos(theta1)*cos(theta2)+(1/2)*sin(theta1)*sin(theta2)*cos(phi2-phi1)

    x = 0
    for l in range(30):
        x += (1j**l)*(2*l+1)*spherical_jn(l,-1*1j*J)*eval_legendre(l,cos(gamma))

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

    theta_list = np.radians(np.linspace(0,180,10)) #Theta(in degrees) (0,180,10)
    phi_list = np.radians(np.linspace(0,342,20)) #Phi(in degrees) (0,342,20)
    #theta_list = np.radians(theta_list)
    angle_set = list(product(theta_list,phi_list)) #Kartezyen Çarpımı

    hamiltonian = []

    for i in range(len(angle_set)):
        for j in range(len(angle_set)):
            
            #hamiltonian.append(J*dot_Product(angle_set[i],angle_set[j]))
            hamiltonian.append(exp_of_dot_Product_sph(angle_set[i],angle_set[j],J))
            
    hamiltonian_max = np.amax(hamiltonian)

    Transfer = []
    num = 0
    for i in range(len(angle_set)):
        row = []
        for j in range(len(angle_set)):
            #row.append(np.exp(hamiltonian[num]-hamiltonian_max))
            row.append(hamiltonian[num]/hamiltonian_max)
            num = num+1
        Transfer.append(row)

    return np.array(Transfer)

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

def flow_to_excel(flow):

    df_list = []
    file_name = 'matrice.xlsx'
    writer = pd.ExcelWriter(file_name)
    for i in range(len(flow)):
        df_list.append(pd.DataFrame(np.real(flow[i])))
        df_list[i].to_excel(writer, sheet_name='RG_NO_{}'.format(i), float_format='%.5f')
    writer.save()

    return True


'''
def main():

    t_list = []
    t_list.append([60,90,120])
    t_list.append([0,90,180])
    t_list.append([45,90,135])
    t_list.append([30,60,120,150])
    t_list.append([15,30,150,165])
    
    for i in range(len(t_list)):

        f = RG_Flow(10,1,t_list[i])
        file_name = 'matrice_{}'.format(i)+'.xlsx'
        flow_to_excel(f, file_name)

    return True
    

        reduced_flow = []
        for i in range(len(Flow_TM)):
            reduced_matrix = []
            for j in range(20):
                row = []
                for k in range(20):
                    row.append(Flow_TM[i][j*10][k*10])
                reduced_matrix.append(row)
            reduced_flow.append(reduced_matrix)
            reduced_flow[i] = np.array(reduced_flow[i]).reshape(20,20)
'''
        
