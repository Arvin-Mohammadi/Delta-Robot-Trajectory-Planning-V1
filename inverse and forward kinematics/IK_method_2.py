# this needs fixing 
# THE JACOBIAN MATRICES
# in this python file i'll attempt to calculate jacobian matrix and sigularities of the Delta Robot

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
from math import sin, cos, acos, asin
from numpy.linalg import inv
import matplotlib.pyplot as plt 
import time

# =================================================================================================
# -- JACOBIAN MATRICES ----------------------------------------------------------------------------
# =================================================================================================

class Jacobian:
	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_triangle_side=0.3464101615, EE_triangle_side=0.2563435195, alpha=[90, 210, 330]):
		
		# initializing the basic geometry and the given data
		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_triangle_side = EE_triangle_side								# the radius of the end effoctor e 
		self.base_triangle_side = base_triangle_side							# the radius of the base or f

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
			t = (-A + (A**2 + B**2 - M**2)**0.5)/(M - B)
			theta_1[i] = 2*math.atan(t)
			theta_2[i] = acos((A - a*cos(theta_1[i]))/(b*sin(theta_3[i]))) - theta_1[i]
		
		theta_1 *= 180/math.pi
		theta_2 *= 180/math.pi
		theta_3 *= 180/math.pi

		theta_1_new = theta_1 - 90 # need this for writing the file 
		
		return (theta_1, theta_2, theta_3)
	
	def get_jacobian_matrix(self, theta_1, thetea_2, theta_3): 
		
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
		py = b - sin(alpha)*(R -r + a*cos(theta_1) + b*sin(theta_3)*cos(theta_2 + theta_1)) + (cos(alpha) - 1/cos(alpha))*b*cos(theta_3) # not correct 
		pz = a*sin(theta_2) + b*sin(theta_3)*sin(theta_2 + theta_1) # not correct 

		return (px, py, pz)
# =================================================================================================
# -- TEST ----------------------------------------------------------------------
# =================================================================================================
t = time.time()

jacobian = Jacobian([0.05, 0.05, -0.31])
print(jacobian.inverse_kinematics())
temp_theta = np.array(jacobian.inverse_kinematics()).transpose()[0]
print(temp_theta)
print(jacobian.forward_kinematics(temp_theta))

jacobian = Jacobian([-0.05, 0.05, -0.40])
print(jacobian.inverse_kinematics())
temp_theta = np.array(jacobian.inverse_kinematics()).transpose()[0]
print(temp_theta)
print(jacobian.forward_kinematics(temp_theta))

jacobian = Jacobian([-0.10, -0.10, -0.42])
print(jacobian.inverse_kinematics())
temp_theta = np.array(jacobian.inverse_kinematics()).transpose()[0]
print(temp_theta)
print(jacobian.forward_kinematics(temp_theta))

jacobian = Jacobian([0.0, -0.15, -0.42])
print(jacobian.inverse_kinematics())
temp_theta = np.array(jacobian.inverse_kinematics()).transpose()[0]
print(temp_theta)
print(jacobian.forward_kinematics(temp_theta))

jacobian = Jacobian([-0.15, 0.1, -0.42])
print(jacobian.inverse_kinematics())
temp_theta = np.array(jacobian.inverse_kinematics()).transpose()[0]
print(temp_theta)
print(jacobian.forward_kinematics(temp_theta))