# THE JACOBIAN MATRICES
# in this python file i'll attempt to calculate jacobian matrix and sigularities of the Delta Robot

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

from re import X
import numpy as np 
import math 
from math import sin, cos 
from numpy.linalg import inv
import matplotlib.pyplot as plt 

# =================================================================================================
# -- JACOBIAN MATRICES ----------------------------------------------------------------------------
# =================================================================================================

# in this section we calculate J_p and J_Theta_dot_dot from the input of: geometric features and position of EE 
# the entire purpose of the jacobian matrix is to relate the Theta_dot_dot matrix to EE velocity matrix 

class Jacobian:
	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):
		# i made a mistake during naming of the variables, base and EE "radius" are not exactly radiuses, they are sides of a triangle ...
		# ... made by EE joints (the joints are the middle of the side of the traingle)
		
		# initializing the basic geometry and the given data
		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_radius = EE_radius								# the radius of the end effoctor e 
		self.base_radius = base_radius							# the radius of the base or f

	def get_theta_ij(self): 
		# this also calculate inverse kinematic

		# assigning constants
		alpha = math.pi/180*self.alpha
		R = self.base_radius*(3**0.5/6)
		r = self.EE_radius*(3**0.5/6)
		a = self.active_rod
		b = self.passive_rod

		# assigning position of the EE 
		px = self.EE_position_global[0]
		py = self.EE_position_global[1]
		pz = self.EE_position_global[2]

		# initializing theta 1, 2, 3
		theta_1 = np.zeros((3))
		theta_2 = np.zeros((3))
		theta_3 = np.zeros((3))

		# calculating theta 1, 2, 3
		for i in [0, 1, 2]:
			theta_3[i] = math.acos((px*sin(alpha[i]) + py*cos(alpha[i]))/b)

			A = px*cos(alpha[i]) - py*sin(alpha[i]) - R + r
			B = pz
			M = (A**2 + B**2 + a**2 - (b*sin(theta_3[i]))**2)/(2*a)
			t = (B + (B**2 - M**2 + A**2)**0.5)/(M + A)
			theta_1[i] = 2*math.atan(t)
			theta_2[i] = math.asin((pz - a*sin(theta_1[i]))/(b*sin(theta_3[i]))) - theta_1[i]
		
		# now we have all of the theta 1, 2 and 3 
		self.theta_1 = theta_1
		self.theta_2 = theta_2
		self.theta_3 = theta_3
		return theta_1
	
	def get_jacobian_matrix(self): 
		
		# initializing J_ij 
		jx = np.zeros((3))
		jy = np.zeros((3))
		jz = np.zeros((3))
		J_P = np.zeros((3, 3))
		J_theta = np.zeros((3, 3))

		for i in [0, 1, 2]:
			jx[i] =  sin(self.theta_3[i])*cos(self.theta_2[i] + self.theta_1[i])*cos(self.alpha[i]) + cos(self.theta_3[i])*sin(self.alpha[i])
			jy[i] = -sin(self.theta_3[i])*cos(self.theta_2[i] + self.theta_1[i])*sin(self.alpha[i]) + cos(self.theta_3[i])*cos(self.alpha[i])
			jz[i] =  sin(self.theta_3[i])*sin(self.theta_2[i] + self.theta_1[i])
			J_P[i, :] = [jx[i], jy[i], jz[i]]
			J_theta[i, i] = sin(self.theta_2[i])*sin(self.theta_3[i])
		
		return (J_P, J_theta)

# =================================================================================================
# -- TEST FOR CONSTANT VELOCITY -------------------------------------------------------------------
# =================================================================================================

# in this part i set up a constant velocity profile as (v = 2 m/s) hense (x = 2*t m)
# and then i discretize it for (n+1) points in space-time

def trapezoidal_generator(p1, p2, n=101, max_velo=0.001): 

	# assigning the first and final EE position 
	p1 = np.array(p1)
	p2 = np.array(p2)

	# calculating a and T 
	T = p1[0]/max_velo*3/2
	a = 3*max_velo/T

	# initilizing matrices 
	position = np.zeros((n+1, 3))
	velocity = np.zeros((n+1, 3))
	vx = np.zeros((n+1))
	vy = np.zeros((n+1))
	vz = np.zeros((n+1))
	x = np.zeros((n+1))
	y = np.zeros((n+1))
	z = np.zeros((n+1))

	# calculating velocity profile 3D
	vx = ???
	vy[:] = 0.0
	vz[:] = 0.0
	velocity[:, 0] = vx
	velocity[:, 1] = vy
	velocity[:, 2] = vz

	# calculating position profile 3D
	x = ???
	y[:] = 0.0
	z[:] = -0.38
	position[:, 0] = x
	position[:, 1] = y
	position[:, 2] = z

	return (position, v)

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

(position, v) = trapezoidal_generator([0.0, 0.0, -0.38], [0.1, 0.0, -0.38])

print("\n\n============== coordinate and velocity ==============\n\n")

print("\n\n this is position matrix \n\n", position, "\n")
print(position.shape, "\n")
print("\n\n this is veloctiy matrix \n\n", v, "\n")
print(v.shape, "\n")

print("\n\n============== jacobian of the velocity ==============\n\n")

Theta_dot = np.zeros((101, 3))

for i in range(0, 101):
	jac = Jacobian(position[i, :])
	jac.get_theta_ij()
	(J_P, J_theta) = jac.get_jacobian_matrix()
	J_theta_inv = inv(J_theta)
	Theta_dot[i, :] = np.matmul(np.matmul(J_theta_inv, J_P), v[i, :])

print("\n\n this is Theta_dot matrix \n\n", Theta_dot)
print(Theta_dot.shape)

plt.title("theta_dot-t plot", fontsize='16')
plt.plot(Theta_dot[:, 0])
plt.plot(Theta_dot[:, 1])
plt.plot(Theta_dot[:, 2])
plt.xlabel("t", fontsize='13')
plt.ylabel("theta_dot", fontsize='13')
plt.legend(('theta_dot_1', 'theta_dot_2', 'theta_dot_3'), loc='best')
plt.grid()
plt.show()
