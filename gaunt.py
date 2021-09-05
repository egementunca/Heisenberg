from sympy.physics.wigner import gaunt
import numpy as np

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

gaunt_coef(10, 5)