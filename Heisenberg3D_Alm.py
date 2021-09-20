import numpy as np
from numpy import cos, sin
from scipy.special import spherical_jn, sph_harm, eval_legendre
from sympy.physics.wigner import gaunt, clebsch_gordan
import pandas as pd
from sympy import re

def renormalize(AlmList):
	flat_list = []
	for i in range(len(AlmList)):
		for j in range(len(AlmList[i])):
			flat_list.append(AlmList[i][j])
	max_val = np.amax(flat_list)
	for i in range(len(AlmList)):
		for j in range(len(AlmList[i])):
			AlmList[i][j] = AlmList[i][j]/max_val
	return AlmList


def lambda_eval(J):
	coef_list = []
	for l in range(11):
		val = (4*np.pi)*(1j**l)*spherical_jn(l,-1j*J)
		coef_list.append(val)

	return np.sum(np.real(coef_list)), np.real(coef_list)

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula
def lambda_eval_formula(J,l):
	val = (4*np.pi)*(1j**l)*spherical_jn(l,-1j*J)
	return val

#Evaluates the Fourier Series Coefficient
def createAlm(J,l,m):
	x = lambda_eval_formula(J, l)*sph_harm(-m, l, 0, np.pi/2)*((-1)**m)
	return x

#Takes a list of lambda values and simply squares for same bonds and multiplies for 2 diff bonds
#Then normalizes it  by dividing list to the biggest element of list

def decimation(AlmList):

	decimatedAlmList = []
	for l in range(len(AlmList)):
		m_list = []
		for m in range(len(AlmList[l])):
			val = (AlmList[l][m])**2
			m_list.append(val)
		decimatedAlmList.append(m_list)

	flat_list = []
	for i in range(len(decimatedAlmList)):
		for j in range(len(decimatedAlmList[i])):
			flat_list.append(decimatedAlmList[i][j])
	max_val = np.amax(flat_list)
	for i in range(len(decimatedAlmList)):
		for j in range(len(decimatedAlmList[i])):
			decimatedAlmList[i][j] = decimatedAlmList[i][j]/max_val

	return decimatedAlmList

def Almprime(J1, J2, l_prec):
	AlmprimeList = []
	for l in range(l_prec):
		m_list =[]
		for m in range(-l,l+1):
			val = 0
			for l1 in range(l_prec):
				for l2 in range(l_prec):
					for m1 in range(-l1, l1+1):
						for m2 in range(-l2, l2+1):
							if (m1+m2) == m:
								if (l1+l2+l)%2 == 0:
									if abs(l1-l2) <= l and l1+l2 >= l:
										x = createAlm(J1, l1, m1)*createAlm(J2, l2, m2)*(((-1)**m))*(gaunt(l1,l2,l,m1,m2,-m).n(8))
										val = val + re(x)
									else:
										pass
								else:
									pass
							else:
								pass
			m_list.append(val)
		AlmprimeList.append(m_list)
	AlmprimeList = renormalize(AlmprimeList)

def AlmBondMove(AlmList1, AlmList2, l_prec):
	AlmBondMoveList = []
	for l in range(l_prec):
		m_list =[]
		for m in range(-l, l+1):
			val = 0
			for l1 in range(l_prec):
				for l2 in range(l_prec):
					for m1 in range(-l1, l1+1):
						for m2 in range(-l2, l2+1):
							if (m1+m2) == m:
								if (l1+l2+l)%2 == 0:
									if abs(l1-l2) <= l and l1+l2 >= l:
										x = AlmList1[l1][m1+l1]*AlmList2[l2][m2+l2]*(((-1)**m))*(gaunt(l1,l2,l,m1,m2,-m).n(8))
										val = val + re(x)
									else:
										pass
								else:
									pass
							else:
								pass
			m_list.append(val)
		AlmBondMoveList.append(m_list)
	AlmBondMoveList = renormalize(AlmBondMoveList)
	return AlmBondMoveList
	
def Renorm_Group_Initial(J, l_prec):

	AlmPrimeList = Almprime(J, J, l_prec)
	AlmDoublePrimeList = AlmBondMove(AlmPrimeList, AlmPrimeList, l_prec)
	AlmTransformed = decimation(AlmDoublePrimeList)
	return AlmTransformed

def Renorm_Group(AlmList1, AlmList2, l_prec):

	AlmPrimeList = AlmBondMove(AlmList1, AlmList2, l_prec)
	AlmDoublePrimeList = AlmBondMove(AlmPrimeList, AlmPrimeList, l_prec)
	AlmTransformed = decimation(AlmDoublePrimeList)
	return AlmTransformed

def flow_to_excel(flow, file_name):

	#Names Index
	row_names = []
	for l in range(len(flow[0])):
		for m in range(-l,1):
			name = "A{}{}".format(l,abs(m))
			row_names.append(name)
	col_names = []
	for n in range(len(flow)-1):
		name = "RG_{}".format(n+1)
		col_names.append(name)
	col_names = ["J = 10"] + col_names

	#Lists Alm through y axis and Renorm group through x axis
	df_list = np.zeros(28*len(flow)).reshape(28,len(flow))
	for x in range(len(flow)):
		flat = []
		for l in range(len(flow[x])):
			for m in range(int((len(flow[x][l])/2)+.5)):
				flat.append(flow[x][l][-m-1])
		for y in range(28):
			df_list[y,x] = flat[y]

	#saving options
	writer = pd.ExcelWriter(file_name)
	df = pd.DataFrame(df_list)
	df.columns = col_names
	df.index = row_names
	df.to_excel(writer, sheet_name='RG', float_format='%1.4f')
	writer.save()
	return True

def flowCreator(J, step_no, l_prec):

	#Alm's created through rg steps stored in flow
	flow = []	

	#q stores original Alm's
	q = []
	for l in range(7):
		m_list = []
		for m in range(-l,l+1):
			val = np.real(createAlm(J, l, m))
			m_list.append(val)
		q.append(m_list)	
	#Renormalizes q
	flat_q = []
	for i in range(len(q)):
		for j in range(len(q[i])):
			flat_q.append(q[i][j])
	max_val = np.amax(flat_q)
	for i in range(len(q)):
		for j in range(len(q[i])):
			q[i][j] = q[i][j]/max_val	
	flow.append(q)	
	
	#First Renorm Group
	a = Renorm_Group_Initial(J, l_prec)
	flow.append(a)
	#Recursive Renorm Group
	for i in range(step_no):
		b = Renorm_Group(a, a, l_prec)
		flow.append(b)
		a = b

	return flow

def threecos(J,theta,phi):

	val = np.exp(J*sin(theta)*cos(phi))
	return val

def expansion(J, theta, phi, l_prec):
	x = 0
	for l in range(l_prec):
		for m in range(-l, l+1):
			x += createAlm(J, l, m)*sph_harm(m,l,phi,theta)
	return x

"""
thetas  = np.linspace(0,np.pi,20)
phis = np.linspace(0,2*np.pi,20)
J = 10
for theta in thetas:
	for phi in phis:
		diff = abs(threecos(J, theta, phi)-expansion(J, theta, phi, ))
		print(diff)"""


flow = flowCreator(10, 10, 25)
flow_to_excel(flow, "lowtemp.xlsx")