import numpy as np
from scipy.special import spherical_jn
from sympy.physics.wigner import gaunt
import pandas as pd

np.random.seed(17)

#since gaunt's integral will be used repeatedly we calculate once and for all for our next calculations
def gaunt_vals(l1,l2):
	val = 0
	for l3 in range(abs(l2-l1),l1+l2+1):
		for m1 in range(-l1,l1+1):
			for m2 in range(-l2,l2+1):
				for m3 in range(-l3,l3+1):
					if (l1+l2+l3)%2 == 0:
						if m3 == m1+m2:
							x = gaunt(l1,l2,l3,m1,m2,-m3,prec=8)*((-1)**m3)
							val += x
						else:
							pass
					else:
						pass
	return val
	
gaunt_table = np.zeros(121).reshape(11,11)
for l1 in range(0,11):
	for l2 in range(0,11):
		gaunt_table[l1,l2] = gaunt_vals(l1, l2)

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula (with sum)
def lambda_eval(J):
	coef_list = []
	for l in range(11):
		val = (2*l+1)*(1j**l)*spherical_jn(l,-1j*J)
		coef_list.append(val)
	return np.sum(coef_list), coef_list

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula  (only formula)
def lambda_eval_formula(J,l):
	val = (2*l+1)*(1j**l)*spherical_jn(l,-1j*J)
	return val

#Takes a list of lambda values and simply squares for same bonds and multiplies for 2 diff bonds
#Then, divides list to the biggest element of list and normalizes it
def decimation(Alm_list):

	decimated_Alm_list = np.zeros(121).reshape(11,11)
	for l1 in range(11):
		for l2 in range(11):
			decimated_Alm_list[l1,l2] = (Alm_list[l1,l2]**2)

	decimated_Alm_list = decimated_Alm_list / np.amax(decimated_Alm_list)

	return decimated_Alm_list

#Evaluates the Coefficient A_lm in Laplace Series for Bond Moving (with sum)
def eval_Alm(J):

	coef_sum, coef_list = lambda_eval(J)
	coef_Alm_list = []
	for l in range(11):
		val = coef_list[l]*(4*np.pi)*(1/(2*l+1))*np.sqrt((2*l+1)*(1/(4*np.pi)))
		coef_Alm_list.append(val)
	return np.sum(coef_Alm_list), coef_Alm_list

#Evaluates the Coefficient A_lm in Laplace Series for Bond Moving (only formula)
def eval_Alm_formula(J,l):
	val = lambda_eval_formula(J, l)*np.sqrt((4*np.pi)*(1/(2*l+1)))
	return val

def eval_iterated_Alm(J1, J2):
	iterated_Alm = np.zeros(121,dtype=np.cdouble).reshape(11,11)
	for l1 in range(0,11):
		for l2 in range(0,11):
			element = eval_Alm_formula(J1, l1)*eval_Alm_formula(J2, l2)*gaunt_table[l1,l2]
			iterated_Alm[l1,l2] = element
	iterated_Alm = np.real(iterated_Alm)
	return iterated_Alm

def bond_move(A_lm1, A_lm2):

	iterated_Alm = np.zeros(121,dtype=np.cdouble).reshape(11,11)
	for l1 in range(0,11):
		for l2 in range(0,11):
			element = np.sum(A_lm1)*np.sum(A_lm2)*gaunt_table[l1,l2]
			iterated_Alm[l1,l2] = element
	iterated_Alm = np.real(iterated_Alm)
	return iterated_Alm

# 2 Bond Moving and 1 Decimation Process for 3 dimensional Renorm Group
def Renorm_Group_Initial(J):

	A_prime = eval_iterated_Alm(J, J)
	A_double_prime = bond_move(A_prime, A_prime)
	A_transformed = decimation(A_double_prime)

	return A_transformed

def Renorm_Group_Transform(Alm_list):

	A_prime = bond_move(Alm_list, Alm_list)
	A_double_prime = bond_move(A_prime, A_prime)
	A_transformed = decimation(A_double_prime)

	return A_transformed

def flow_to_excel(flow, file_name):

    df_list = []
    writer = pd.ExcelWriter(file_name)
    for i in range(len(flow)):
        df_list.append(pd.DataFrame(flow[i]))
        df_list[i].to_excel(writer, sheet_name='RG_NO_{}'.format(i+1), float_format='%1.5f')
    writer.save()

    return True


flow = []
a = Renorm_Group_Initial(.01)
flow.append(a)
for i in range(5):
	b = Renorm_Group_Transform(a)
	flow.append(b)
	print(b)
	a = b





"""
array = np.real(eval_iterated_Alm(10, 10))
np.savetxt("foo_low.csv",array,delimiter=',',fmt="%1.5e")
"""










