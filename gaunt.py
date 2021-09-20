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

def gaunt_coef(l,m):

	vals = []
	for l1 in range(6):
		for l2 in range(6):
			for m1 in range(-l1,l1+1):
				for m2 in range(-l2,l2+1):
					x = gaunt(l1,l2,l,m1,m2,-m,prec=8)*((-1)**m)
					print(x)
					vals.append(x)
	return vals

def flow_to_excel(flow, file_name):

    df_list = []
    writer = pd.ExcelWriter(file_name)
    for i in range(len(flow)):
    	df_list.append(flow[i])
    df = pd.DataFrame(df_list)
    df.to_excel(writer, sheet_name='RG', float_format='%1.4f')
    writer.save()

    return True

flow = []
a = gaunt_coef(1, -1)
print(sum(a))
#for l in range(6):
#	for m in range(-l,l+1):
#		flow.append(gaunt_coef(l, m))
#flow_to_excel(flow, "gaunt.xlsx")

