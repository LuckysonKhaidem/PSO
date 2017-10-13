import pso
import numpy as np
def f1(x):
	x = np.array(x)
	return (x**2).sum()

sol = pso.optimize(objective_function = f1,swarm_size = 500,tolerance = 0.005,number_of_dimensions = 10,omega = 0.05,c1 = 0.05,c2 = 0.85,max_iteration =1000)
print "Solution is", sol
