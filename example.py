import pso
import numpy as np
import math

# Minimize f(x) = x^2 + y^2
# where -inf < x < inf,
#	-inf < y < inf
# such that x + y >= 10


def penalty(k,x):
	if x > 0:
		return k*x**2
	else:
		return 0
def f1(x):
	return (x[0]**2 + x[1]**2) + penalty(100,10 - (x[0] + x[1]) )


sol = pso.optimize(objective_function = f1,swarm_size = 5000,tolerance = 0.005,number_of_dimensions = 2,omega = 0.05,c1 = 0.05,c2 = 0.85,lower_bound = [-99999,-99999], upper_bound = [99999,99999],max_iteration =1000)
print "Solution is", sol
