from sympy.physics.wigner import gaunt
x=0

for l3 in range(0,10):
	for l2 in range(0,10):
		for l1 in range(0,10):
			for m2 in range(-l1,l1+1):
				for m1 in range(-l2,l2+1):
					if ((l1+l2+l3)%2) == 0:
						if (abs(l1-l2) <= l3) and (l3 <= l1+l2):
							if (abs(-m1-m2) <= l3):
								val = gaunt(l1,l2,l3,m1,m2,-m1-m2)
								x = x + val
								print(val)
							else:
								pass
						else:
							pass
					else:
						pass
print(x.n(20))