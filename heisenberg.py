import random
import numpy as np
from numpy import sin, cos
from itertools import product

# Dot product of unit vector spins s_1 * s_2
def dot_Product(angle_set1, angle_set2):

    theta1, phi1 = angle_set1[0], angle_set1[1]
    theta2, phi2 = angle_set2[0], angle_set2[1]
    dot_product = sin(theta1)*sin(theta2)*cos(phi1)*cos(phi2)+sin(theta1)*sin(theta2)*sin(phi1)*sin(phi2)+cos(theta1)*cos(theta2)

    return dot_product

#Probably wrong!!!, an estimation for Transfer Matrix for Heisenberg Model
#exp(i*s_i*s_j) = (Spherical Harmonics)
def TransferMatrixCreator(J):

    theta_list = np.radians(np.linspace(0,180,7)) #Theta(in degrees) : [0, 30, 60, 90, 120, 150, 180]
    phi_list = np.radians(np.linspace(0,330,12)) #Phi(in degrees) : [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    angle_set = list(product(theta_list,phi_list)) #Kartezyen Çarpımı ile 84 elemanlı olabilecek theta, phi kombinasyonları

    hamiltonian = []

    for i in range(len(angle_set)):
        for j in range(len(angle_set)):
            hamiltonian.append(J*dot_Product(angle_set[i],angle_set[j]))

    Transfer = []
    num = 0
    for i in range(len(angle_set)):
        row = []
        for j in range(len(angle_set)):
            row.append(np.exp(hamiltonian[num]))
            num = num+1
        Transfer.append(row)

    return np.array(Transfer)

def TransferMatrixDec(matrix1, matrix2):

    Transfer = np.dot(matrix1,matrix2)
    return Transfer

def RenormalizationGroup(bondNumber, transfer_matrix_list):
    
    np.random.seed(17)
    transformed_matrix_list = []

    for i in range(bondNumber):

        Transfer = []
        Transfer = TransferMatrixDec(transfer_matrix_list[np.random.randint(1,bondNumber)],transfer_matrix_list[np.random.randint(1,bondNumber)])
        transformed_matrix_list.append(Transfer)

    return transformed_matrix_list

def RG_Flow(bondNumber, RG_step, J_initial):

    bN = bondNumber
    tm_list = []

    for i in range(bN):
        tm = TransferMatrixCreator(J_initial)
        tm_list.append(tm)

    for j in range(RG_step):

        tm_list_transformed = RenormalizationGroup(bN, np.array(tm_list))

#################SHOULD CHECK SMTHNG   


#################TO UNDRSTND CALCULATIONS
        
        tm_list = tm_list_transformed

        print(tm_list)

    return True

RG_Flow(500,20,1)


        











