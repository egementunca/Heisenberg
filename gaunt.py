from sympy.physics.wigner import gaunt, clebsch_gordan
import numpy as np
import pandas as pd

def cleb_coef(l):
	vals = []
	for l1 in range(1000):
		for l2 in range(1000):
			if l <= l1+l2 and l >= abs(l1-l2):
				if (l1+l2+l)%2 == 0:
					x = (clebsch_gordan(l1,l2,l,0,0,0).n(8))**2
					vals.append(x)
					print("l = {}, l1 = {}, l2 = {} value is {}".format(l,l1,l2,x))
				else:
					pass
			else:
				pass
	return(np.sum(vals))

a = cleb_coef(1)
print(a)

def gaunt_coef(l1,l2):

	vals = []
	count = 0
	for l3 in range(abs(l2-l1),l1+l2+1):
		for m1 in range(-l1,l1+1):
			for m2 in range(-l2,l2+1):
				for m3 in range(-l3,l3+1):
					if m1+m2 == m3:
						if (l1+l2+l3)%2 == 0:
							count = count+1
							x = gaunt(l1,l2,l3,m1,m2,-m3,prec=8)*((-1)**m3)
							print("Gaunt Integral for l3={}, m1={}, m2={}, m3={} is {}".format(l3,m1,m2,m3,x))
							vals.append(x)
						else:
							#print("Gaunt Integral for l3={}, m1={}, m2={}, m3={} cannot be evaluated because of odd L".format(l3,m1,m2,m3))
							pass
					else:
						#print("Gaunt Integral for l3={}, m1={}, m2={}, m3={} cannot be evaluated because of non zero M".format(l3,m1,m2,m3))
						pass
	print("Gaunt Integral for l1={} and l2={} summed over {} possible L,Ms".format(l1,l2,count))
	print(np.sum(vals))
	return np.sum(vals)

def flow_to_excel(flow, file_name):

    df_list = []
    writer = pd.ExcelWriter(file_name)
#len(flow)
    for i in range(1):
        df_list.append(pd.DataFrame(flow))
        df_list[i].to_excel(writer, sheet_name='RG_NO_{}'.format(i+1), float_format='%1.5f')
    writer.save()

    return True

