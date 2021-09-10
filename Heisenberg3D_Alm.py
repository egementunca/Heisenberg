import numpy as np
from scipy.special import spherical_jn, sph_harm
from sympy.physics.wigner import gaunt, clebsch_gordan
import pandas as pd
from sympy import re

np.random.seed(17)

def lambda_eval(J):
	coef_list = []
	for l in range(11):
		val = (4*np.pi)*(1j**l)*spherical_jn(l,-1j*J)
		coef_list.append(val)
	return np.sum(np.real(coef_list)), np.real(coef_list)

#Evaluates the Coefficient Lambda in Plane Wave Expansion Formula  (only formula)
def lambda_eval_formula(J,l):
	val = (4*np.pi)*(1j**l)*spherical_jn(l,-1j*J)
	return val

def createAlm(J,l,m):
	x = lambda_eval_formula(J, l)*sph_harm(m,l,0,0)
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

def Almprime(J1, J2):
	AlmprimeList = []
	for l in range(6):
		m_list =[]
		for m in range(-l,l+1):
			val = 0
			for l1 in range(6):
				for l2 in range(6):
					for m1 in range(-l1, l1+1):
						for m2 in range(-l2, l2+1):
							if abs(l1-l2) <= l and l1+l2 >= l:
								if (l1+l2+l)%2 == 0:
									if (m1+m2) == m:
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
	flat_list = []
	for i in range(len(AlmprimeList)):
		for j in range(len(AlmprimeList[i])):
			flat_list.append(AlmprimeList[i][j])
	max_val = np.amax(flat_list)
	for i in range(len(AlmprimeList)):
		for j in range(len(AlmprimeList[i])):
			AlmprimeList[i][j] = AlmprimeList[i][j]/max_val

	return AlmprimeList

def AlmBondMove(AlmList1, AlmList2):
	AlmBondMoveList = []
	for l in range(6):
		m_list =[]
		for m in range(-l, l+1):
			val = 0
			for l1 in range(6):
				for l2 in range(6):
					for m1 in range(-l1, l1+1):
						for m2 in range(-l2, l2+1):
							if abs(l1-l2) <= l and l1+l2 >= l:
								if (l1+l2+l)%2 == 0:
									if (m1+m2) == m:
										x = AlmList1[l1][m1+l1]*AlmList2[l2][m2+l2]*(((-1)**m))*gaunt(l1,l2,l,m1,m2,-m).n(8)
										val = val + re(x)
									else:
										pass
								else:
									pass
							else:
								pass
			m_list.append(val)
		AlmBondMoveList.append(m_list)
	flat_list = []
	for i in range(len(AlmBondMoveList)):
		for j in range(len(AlmBondMoveList[i])):
			flat_list.append(AlmBondMoveList[i][j])
	max_val = np.amax(flat_list)
	for i in range(len(AlmBondMoveList)):
		for j in range(len(AlmBondMoveList[i])):
			AlmBondMoveList[i][j] = AlmBondMoveList[i][j]/max_val

	return AlmBondMoveList
	

def Renorm_Group_Initial(J):

	AlmPrimeList = Almprime(J, J)
	AlmDoublePrimeList = AlmBondMove(AlmPrimeList,AlmPrimeList)
	AlmTransformed = decimation(AlmDoublePrimeList)
	return AlmTransformed

def Renorm_Group(AlmList1, AlmList2):

	AlmPrimeList = AlmBondMove(AlmList1, AlmList2)
	AlmDoublePrimeList = AlmBondMove(AlmPrimeList, AlmPrimeList)
	AlmTransformed = decimation(AlmDoublePrimeList)
	return AlmTransformed

def flow_to_excel(flow, file_name):

    df_list = []
    writer = pd.ExcelWriter(file_name)
    for x in range(len(flow)):
    	iterated_vals = []
    	for y in range(len(flow[0])):
    		desired_vals = flow[y][x]
    		iterated_vals.append(desired_vals[0:int((len(desired_vals)/2)+.5)])
    	df_list.append(iterated_vals)
    df_list = pd.DataFrame(df_list)
    df_list.to_excel(writer, sheet_name='RG', float_format='%1.5f')
    writer.save()
    return True

def flowCreator(J, step_no):

	#Alm's created through rg steps stored in flow
	flow = []	

	#q stores original Alm's
	q = []
	for l in range(6):
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
	a = Renorm_Group_Initial(J)
	flow.append(a)
	#Recursive Renorm Group
	for i in range(step_no):
		b = Renorm_Group(a,a)
		flow.append(b)
		a = b

	return flow

flow = flowCreator(.1, 4)
flow_to_excel(flow, "rgHeisenberghigh.xlsx")