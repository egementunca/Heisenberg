import numpy as np
import matplotlib.pyplot as plt

np.random.seed(17)

def g(J):

    return 4*np.pi*np.sinh(J)*(1/(J))

def bond_pool_creator(J,size):

	pool = []

	for i in range(size):
		a = np.random.random()
		if a > .5:
			j_val = np.abs(J)*1
		elif a <= .5:
			j_val = np.abs(J)*1
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

	return(lambda_vals)

def main(J, size):

	part_func_iteration = []
	g_vals_iteration = []
	f_en_iteration = []

	bond_pool = bond_pool_creator(J, size)
	part_func = part_func_initial(bond_pool)
	g_vals = np.log(part_func)
	free_energy = np.sum(g_vals)/len(g_vals)
	part_func_iteration.append(part_func)
	g_vals_iteration.append(g_vals)
	f_en_iteration.append(free_energy)

	for i in range(10):

		g_vals_iterated = decimation(g_vals)
		part_func_iterated = np.exp(g_vals_iterated)
		free_energy_new = np.sum(g_vals_iterated)/len(g_vals_iterated)
		part_func_iteration.append(part_func_iterated)
		g_vals_iteration.append(g_vals_iterated)
		f_en_iteration.append(free_energy_new)
		part_func = part_func_iterated
		g_vals = g_vals_iterated

	return part_func_iteration, g_vals_iteration, f_en_iteration

a,b,c = main(1, 1024)
for i in range(11):
	print("part_func len: {}, g_vals_len:{}, free_en:{}").format(len(a[i]),len(b[i]),c[i])
	