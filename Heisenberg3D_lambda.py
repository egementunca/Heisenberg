import numpy as np
from scipy.special import spherical_jn
from sympy.physics.wigner import gaunt, clebsch_gordan
import pandas as pd
from sympy import re

np.random.seed(17)

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula (with sum)
def lambda_eval(J):
	coef_list = []
	for l in range(11):
		val = (2*l+1)*(1j**l)*spherical_jn(l,-1j*J)
		coef_list.append(val)
	return np.sum(np.real(coef_list)), np.real(coef_list)

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula  (only formula)
def lambda_eval_formula(J,l):
	val = (2*l+1)*(1j**l)*spherical_jn(l,-1j*J)
	return val

#Takes a list of lambda values and simply squares for same bonds and multiplies for 2 diff bonds
#Then normalizes it  by dividing list to the biggest element of list
def decimation(lambda_list):

	decimated_lambda_list = []
	for l in range(len(lambda_list)):
		val = ((4*np.pi)/(2*l+1))*(lambda_list[l]**2)
		decimated_lambda_list.append(val)
	decimated_lambda_list = np.array(decimated_lambda_list)/np.amax(decimated_lambda_list)

	return decimated_lambda_list

def lambda_prime(J1, J2):
	lambda_prime_list = []
	for l in range(11):
		val = 0
		for l1 in range(11):
			for l2 in range(11):
				if abs(l1-l2) <= l and l1+l2 >= l:
					if (l1+l2+l)%2 == 0:
						x = lambda_eval_formula(J1, l1)*lambda_eval_formula(J2, l2)*(clebsch_gordan(l1,l2,l,0,0,0)**2)
						val = val + re(x)
					else:
						pass
				else:
					pass
		lambda_prime_list.append(val)
	lambda_prime_list = np.real(lambda_prime_list)
	return lambda_prime_list/np.amax(lambda_prime_list)

def lambda_bond_move(lambda1, lambda2):
	lambda_prime_list = []
	for l in range(11):
		val = 0
		for l1 in range(11):
			for l2 in range(11):
				if abs(l1-l2) <= l and l1+l2 >= l:
					if (l1+l2+l)%2 == 0:
						x = lambda1[l1]*lambda2[l2]*(clebsch_gordan(l1,l2,l,0,0,0)**2)
						val = val + x
					else:
						pass
				else:
					pass
		lambda_prime_list.append(val)
	lambda_prime_list = np.real(lambda_prime_list)
	return lambda_prime_list/np.amax(lambda_prime_list)

# 2 Bond Moving and 1 Decimation Process for 3 dimensional Renorm Group
def Renorm_Group_Initial(J):

	lambda_prime_list = lambda_prime(J, J)
	lambda_double_prime = lambda_bond_move(lambda_prime_list,lambda_prime_list)
	lambda_transformed = decimation(lambda_double_prime)
	return lambda_transformed

def Renorm_Group_Transform(lambda_list):

	lambda_prime = lambda_bond_move(lambda_list, lambda_list)
	lambda_double_prime = lambda_bond_move(lambda_prime,lambda_prime)
	lambda_transformed = decimation(lambda_double_prime)
	return lambda_transformed

def flow_to_excel(flow, file_name):

	df_list = np.zeros(len(flow)*len(flow[0])).reshape(len(flow[0]),len(flow))
	writer = pd.ExcelWriter(file_name)
	for x in range(len(flow)):
		for y in range(len(flow[x])):
			df_list[y,x] = flow[x][y]
	df = pd.DataFrame(df_list)
	df.to_excel(writer, sheet_name='RG', float_format='%1.4f')
	writer.save()
	return True


def flow_creator(J, step_no):

	flow = []

	q = []
	for l in range(11):
		q.append(np.real(lambda_eval_formula(J, l)))
	q_max = np.amax(q)
	q = np.array(q)
	q = q/q_max
	flow.append(q)

	a = Renorm_Group_Initial(J)
	#a = np.ones(11)
	flow.append(a)
	for i in range(step_no):
		b = Renorm_Group_Transform(a)
		flow.append(b)
		a = b
	
	return flow

flow = flow_creator(10, 4)
flow_to_excel(flow, "testlow.xlsx")

