# this needs fixing 
# THE JACOBIAN MATRICES
# in this python file i'll attempt to calculate jacobian matrix and sigularities of the Delta Robot

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
from math import sin, cos, acos, asin, atan
from numpy.linalg import inv
import matplotlib.pyplot as plt 
from numpy import matmul
import time

# =================================================================================================
# -- JACOBIAN MATRICES ----------------------------------------------------------------------------
# =================================================================================================

class Jacobian:
	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_triangle_side=0.3464101615, EE_triangle_side=0.2563435195, alpha=[0, 120, 240]):
		
		# initializing the basic geometry and the given data
		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_triangle_side = EE_triangle_side				# the triangle side of the end effoctor e 
		self.base_triangle_side = base_triangle_side			# the triangle side of the base or f

	def inverse_kinematics(self): 
		# this also calculate inverse kinematic

		# assigning constants
		alpha = math.pi/180*self.alpha
		R = self.base_triangle_side*(3**0.5/6)
		r = self.EE_triangle_side*(3**0.5/6)
		a = self.active_rod
		b = self.passive_rod

		# assigning position of the EE 
		px = self.EE_position_global[0]
		py = self.EE_position_global[1]
		pz = self.EE_position_global[2]

		# px = -abs(px)
		# py = abs(py)

		# initializing theta 1, 2, 3
		theta_1 = np.zeros((3))
		theta_2 = np.zeros((3))
		theta_3 = np.zeros((3))

		# calculating theta 1, 2, 3
		for i in [0, 1, 2]:
			theta_3[i] = acos((px*sin(alpha[i]) + py*cos(alpha[i]))/b)

			A = px*cos(alpha[i]) - py*sin(alpha[i]) - R + r
			B = pz
			M = -(A**2 + B**2 + a**2 - (b*sin(theta_3[i]))**2)/(2*a)
			t = (-B + (A**2 + B**2 - M**2)**0.5)/(M - A)
			theta_1[i] = 2*atan(t)
			theta_2[i] = asin((pz - a*sin(theta_1[i]))/(b*sin(theta_3[i]))) - theta_1[i]

		
		theta_1 *= 180/math.pi
		theta_2 *= 180/math.pi
		theta_3 *= 180/math.pi

		theta_1_new = np.zeros((3))
		theta_1_new[0] = theta_1[0]  
		theta_1_new[1] = theta_1[1] 
		theta_1_new[2] = theta_1[2] 

		return (theta_1, theta_2, theta_3)
	
	def get_jacobian_matrix(self, theta_1, theta_2, theta_3): 
		
		# changing the angles back to rad
		theta_1 /= 180/math.pi
		theta_2 /= 180/math.pi
		theta_3 /= 180/math.pi

		# initializing J_ij 
		jx = np.zeros((3))
		jy = np.zeros((3))
		jz = np.zeros((3))
		J_P = np.zeros((3, 3))
		J_theta = np.zeros((3, 3))

		for i in [0, 1, 2]:
			jx[i] =  sin(theta_3[i])*cos(theta_2[i] + theta_1[i])*cos(self.alpha[i]) + cos(theta_3[i])*sin(self.alpha[i])
			jy[i] = -sin(theta_3[i])*cos(theta_2[i] + theta_1[i])*sin(self.alpha[i]) + cos(theta_3[i])*cos(self.alpha[i])
			jz[i] =  sin(theta_3[i])*sin(theta_2[i] + theta_1[i])
			J_P[i, :] = [jx[i], jy[i], jz[i]]
			J_theta[i, i] = sin(theta_2[i])*sin(theta_3[i])
		
		return (J_P, J_theta)

	def forward_kinematics(self, theta): # theta is theta_1, theta_2, theta_3 from one of the motor number one

		# assigning constants
		alpha = math.pi/180*self.alpha[0]
		R = self.base_triangle_side*(3**0.5/6)
		r = self.EE_triangle_side*(3**0.5/6)
		a = self.active_rod
		b = self.passive_rod
		theta *= math.pi/180
		theta_1 = theta[0]
		theta_2 = theta[1]
		theta_3 = theta[2]

		px = cos(alpha)*(R - r + a*cos(theta_1) + b*sin(theta_3)*cos(theta_2 + theta_1)) + sin(alpha)*b*cos(theta_3)
		py = sin(alpha)*(r - R - a*cos(theta_1) - b*sin(theta_3)*cos(theta_2 + theta_1)) + cos(alpha)*b*cos(theta_3)
		pz = a*sin(theta_1) + b*sin(theta_3)*sin(theta_2 + theta_1) # not correct why ???
		
		return (px, py, pz)

	def theta_dot(self, j_p, j_theta, v):
		return matmul(matmul(inv(j_theta), j_p), v)


# =================================================================================================
# -- TEST ----------------------------------------------------------------------
# =================================================================================================

t = np.linspace(0, 1, 1000)
gamma = t*2*np.pi 
x = 0.1*np.cos(gamma)
y = 0.1*np.sin(gamma) 
z = -0.35

x_dot = 0.1*2*np.pi*(-np.sin(gamma))
y_dot = 0.1*2*np.pi*( np.cos(gamma))


THETA_DOT = np.zeros((1000, 3))
for i, t_i in enumerate(t):
	jacobian = Jacobian([x[i], y[i], z])
	theta_1, theta_2, theta_3 = jacobian.inverse_kinematics()
	j_p, j_theta = jacobian.get_jacobian_matrix(theta_1, theta_2, theta_3)
	THETA_DOT[i, :] = jacobian.theta_dot(j_p, j_theta, [x_dot[i], y_dot[i], 0])

plt.plot(t, THETA_DOT[:, :])
plt.show()
