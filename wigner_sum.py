from sympy.physics.wigner import gaunt
import numpy as np
x=0

for l3 in range(0,4):
	for l2 in range(0,4):
		for l1 in range(0,4):
			for m2 in range(-l2,l2+1):
				for m1 in range(-l1,l1+1):
					if ((l1+l2+l3)%2) == 0:
						if (abs(l1-l2) <= l3) and (l3 <= l1+l2):
							if (abs(-m1-m2) <= l3):
								val = np.sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)*(1/(4*np.pi)))*(-1**(m1+m2))*gaunt(l1,l2,l3,m1,m2,-m1-m2).n(10)
								x = x + val
								print("integral of l1={},l2={},l3={},m1={},m2={},m3={} is {}".format(l1,l2,l3,m1,m2,(-m1-m2),val))
							else:
								print("cant compute the integral of l1={},l2={},l3={},m1={},m2={},m3={}".format(l1,l2,l3,m1,m2,(-m1-m2)))
						else:
							print("cant compute the integral of l1={},l2={},l3={},m1={},m2={},m3={}".format(l1,l2,l3,m1,m2,(-m1-m2)))
					else:
						print("cant compute the integral of l1={},l2={},l3={},m1={},m2={},m3={}".format(l1,l2,l3,m1,m2,(-m1-m2)))
print(x)