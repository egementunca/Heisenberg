import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, eval_legendre, spherical_in

np.random.seed(17)

rndm_list = []
for i in range(1024):
	rndm_list.append(np.random.random())

def g(J):

    return 4*np.pi*np.sinh(J)*(1/(J))

def V_tilde(J):

	return np.exp(-J) * spherical_in(J)

def coefficient(l,J):
	return 4*np.pi*(1j**l)*spherical_jn(l,-1j*J)


def V_tilde_prime(J, s):
	x = 0
	for p in range(-10,11):
		try:
			val = V_tilde(J, p)*V_tilde(J, s-p)
		except:
			val =0
		x += val
	return x**2

def V_prime(J, theta):
	x = 0
	for s in range(-10, 11):
		val = (np.exp(-1j*s*theta))*V_tilde_prime(J, s)
		x += val
	return x

print(V_prime(1, 0))

def bond_pool_creator(J, p, size):

	pool = []

	for i in range(size):
		a = rndm_list[i]
		if a > p:
			j_val = np.abs(J)*0.5
		elif a <= p:
			j_val = np.abs(J)*1.5
		pool.append(j_val)

	return pool

def part_func_initial(bond_pool):

	lambda_vals = []
	for i in range(len(bond_pool)):
		lambda_val = g(bond_pool[i])
		lambda_vals.append(lambda_val)
	return lambda_vals
  
def decimation(g_vals):

	lambda_vals = []
	for i in range(len(g_vals)/2):
		lambda_iterated = g_vals[i]+g_vals[i+1]
		lambda_vals.append(lambda_iterated)

	return lambda_vals

def bond_moving():
	
	pass

def main(J, p, size):

	g_vals_iteration = []
	f_en_iteration = []

	bond_pool = bond_pool_creator(J, p, size)
	part_func = part_func_initial(bond_pool)
	g_vals = np.log(np.array(part_func))
	free_energy = np.sum(g_vals)/len(g_vals)
	g_vals_iteration.append(g_vals)
	f_en_iteration.append(free_energy)
	for i in range(10):
		g_vals_iterated = decimation(g_vals)
		free_energy_new = np.sum(g_vals_iterated)/(len(g_vals_iterated))
		g_vals_iteration.append(g_vals_iterated)
		f_en_iteration.append(free_energy_new)
		g_vals = g_vals_iterated

	return g_vals_iteration, f_en_iteration

def cv_data(p):

	eps = 1e-5

	x = np.linspace(0.1,10,500)
	x_right = x+eps
	x_left = x-eps
	y, y_right, y_left = [], [], []
	cv = []
	real_x = []
	for i in range(len(x)):
		val = x[i]*p*1.5+x[i]*(1-p)*0.5
		real_x.append(val)
		y.append(main(x[i],p,1024))
		y_right.append(main(x_right[i],p,1024))
		y_left.append(main(x_left[i],p,1024))
	for i in range(len(x)):
		cv.append((real_x[i]**2)*(y_right[i]+y_left[i]-2*y[i])/(eps**2))


	return real_x, cv
"""
labels = []
p_vals = np.arange(0,11)/10
cv_list = []
x_list = []
for i in range(len(p_vals)):
	x, cv = cv_data(p_vals[i])
	cv_list.append(cv)
	x_list.append(1/np.array(x))
	labels.append("p1:{}, p2:{}".format(p_vals[i],np.around(1-p_vals[i],decimals=1)))
plt.plot(x_list[0],cv_list[0],
		x_list[1],cv_list[1],
		x_list[2],cv_list[2],
		x_list[3],cv_list[3],
		x_list[4],cv_list[4],
		x_list[5],cv_list[5],
		x_list[6],cv_list[6],
		x_list[7],cv_list[7],
		x_list[8],cv_list[8],
		x_list[9],cv_list[9],
        x_list[10],cv_list[10])
plt.legend(labels)
plt.show()
"""