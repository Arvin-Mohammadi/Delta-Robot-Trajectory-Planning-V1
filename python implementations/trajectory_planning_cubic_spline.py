# in this file i attempt making a trajectory generator which makes the Delta robot go around a circle 
# the boundary conditions are as followed: 
# 1. continuity on the position 	--> 2n 		conditions
# 2. continuity on the velocity 	--> n - 1 	conditions 
# 3. continuity on the acceleration --> n - 1 	conditions 
# 4. intial and final velocity zero --> 2 		conditions

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

from http.cookiejar import LWPCookieJar
import time
from typing import final
t = time.time()

import numpy as np 
import math 
import matplotlib.pyplot as plt 

print("time is:")
print(time.time() - t)


# =================================================================================================
# -- FORWARD KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================

# here we have the angles of the joints (theta1, theta2, theta3)
# our goal is to find the position of the end effector (EE center position)

class ForwardKinematics:

	def __init__ (self, theta):
		self.theta = theta 


# =================================================================================================
# -- INVERSE KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================

# here we have the EE position and the lengths of rods and basic geometry
# our goal is to find theta1, theta2 and theta3 

class InverseKinematics:

	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):
		# initializing the basic geometry and the given data

		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_radius = EE_radius								# the radius of the end effoctor e 
		self.base_radius = base_radius							# the radius of the base or f
		

	def get_J1_positions(self):
		
		# initializing the J1, F1 position
		self.F1_position = ([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
		self.J1_position = ([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

		for i in [0, 1, 2]:

			alpha = self.alpha[i]

			# converting the EE position to local position based on alpha
			x = float(self.EE_position_global[0])*math.cos(alpha*math.pi/180) + float(self.EE_position_global[1])*math.sin(alpha*math.pi/180)
			y = - float(self.EE_position_global[0])*math.sin(alpha*math.pi/180) + float(self.EE_position_global[1])*math.cos(alpha*math.pi/180)
			z = float(self.EE_position_global[2])

			self.EE_position_local = np.array([x, y, z])

			# here positions of the the 3 important points are found (E1, E1_prime, F1) localy 
			tan_30_deg = 0.5773502692
			self.E1_position = self.EE_position_local + np.array([0, -self.EE_radius/2*tan_30_deg, 0]) 
			self.E1_prime_position = np.array([0, self.E1_position[1], self.E1_position[2]])
			self.F1_position[i] = np.array([0, - self.base_radius/2*tan_30_deg, 0])

			# and from those 3 point we find J1 
			# (intersection of 2 circles, with centers of F1 and E1_prime and ridus of r_f and (r_e**2 - E1_x**2)**0.5)
			rf = float(self.active_rod)
			re = float(self.passive_rod)
			x0 = float(self.E1_position[0])
			y0 = float(self.E1_position[1])
			z0 = float(self.E1_position[2])
			yF = float(self.F1_position[i][1])

			c1 = (x0**2 + y0**2 + z0**2 + rf**2 - re**2 - yF**2)/(2*z0)
			c2 = (yF - y0)/z0
			c3 = -(c1 + c2*yF)*(c1 + c2*yF) + rf*(c2**2*rf + rf)
			if c3 < 0:
				print("non existing point")
				self.J1_position = 0
				return

			y = (yF - c1*c2 - c3**0.5)/(c2**2 + 1)
			z = c1 + c2*y

			self.J1_position[i] = [0, y, z]
	
	def get_theta(self):
		self.theta = np.zeros((3))
		for i in [0, 1, 2]:
			z_J1 = self.J1_position[i][2]
			y_J1 = self.J1_position[i][1]
			y_F1 = self.F1_position[i][1]

			self.theta[i] = (math.atan(-z_J1/(y_F1 - y_J1)))

		# Gearbox ratio = 50, zero offset = 47.2 degree
		self.theta = ((np.array(self.theta)*180/math.pi) + 47.2)*50

		return self.theta

# =================================================================================================
# -- POSITION GENERATOR ---------------------------------------------------------------------------
# =================================================================================================

# in this part we try to generate a circle and get n sample points from that circle.

class PositionGenerator:

	def __init__(self, ratio, center, n=100):
		pi = math.pi
		self.n = n
		self.ratio = ratio		# r
		self.center = center	# xc, yc, zc
		self.gamma = np.linspace(0, 2*pi, num=self.n)

	def cart_position(self):
		# self.points is the positions of the n points in x, y and z directions (n*3 matrix)
		self.points = np.array([np.cos(self.gamma)*self.ratio + np.ones((self.n))*self.center[0], np.sin(self.gamma)*self.ratio + np.ones((self.n))*self.center[1], np.ones((self.n))*self.center[2]])
		self.points = np.transpose(self.points)
		return self.points
# =================================================================================================
# -- POLYNOMIAL COEFFICIENT MATRIX ----------------------------------------------------------------
# =================================================================================================

# using the boundary conditions we find a[k][i] for k = 0, ..., n and i = 0, 1, 2, 3

class Coeff:

	def __init__(self, points):

		# initialzing the postion, velocity and acceleration vectors 
		# all of the vectors are n*3 matrices (e.g: for v we have 3 components in x, y and z direction)
		self.points = np.array(points)
		self.velo = np.zeros(self.points.shape)
		self.t = np.zeros(self.points.shape)
		self.n = self.points.shape[0]

		# time value assignment, shape = n*3
		for i in [0, 1, 2]:
			self.t[:, i] = np.linspace(0, self.n, num=self.n)
		
		# time intervals (T), shape = n-1*3
		self.T = self.t[1:self.n, :] - self.t[0:self.n-1, :]

	def velocity(self, initial_velo=[0, 0, 0], final_velo=[0, 0, 0]):
		# initializing c_prime and A_prime matrices (last dimension is for the x, y, z directions)
		A_prime = np.zeros((3, self.n-2, self.n-2))
		c_prime = np.zeros((self.n-2, 3))
		T = self.T
		
		# makin the A_prime and c_prime matrices 
		for i in range(A_prime.shape[1]):
			if i != 0:
				A_prime[:, i, i-1] = T[i+1]
			A_prime[:, i, i] = 2*(T[i] + T[i+1])
			if i != A_prime.shape[1]-1:
				A_prime[:, i, i+1] = T[i]

		for i in range(A_prime.shape[1]):
			c_prime[i, :] = 3/(T[i]*T[i+1])*(T[i]**2*(self.points[i+2, :] - self.points[i+1, :]) + T[i+1]**2*(self.points[i+1, :] - self.points[i, :]))
			if i == 0:
				c_prime[i, :] -= T[i+1]*initial_velo
			elif i == A_prime.shape[1]-1:
				c_prime[i, :] -= T[i]*final_velo

		# calculating v vector from A_prime and C_prime matrices 
		v = np.zeros((self.n-2, 3))
		for i in [0, 1, 2]:
			M = np.linalg.inv(A_prime[i, :, :])
			N = c_prime[:, i]
			v[:, i] = np.matmul(M, N)

		self.velo[0, :] = initial_velo
		self.velo[self.n-1, :] = final_velo
		self.velo[1:self.n-1, :] = v

		print("this is A prime\n", A_prime)
		print("this is c prime\n", c_prime)
		print("this is velocity matrix\n", self.velo)

	def coeff_matrix(self):
		# initializing the coefficient matrix 
		# dim 1 == number of polynomials 				--> k = 0, ..., n-2
		# dim 2 == number of x, y, z directions 		--> 3
		# dim 3 == number of coeff in the polynomial 	--> (e.g: a[k][0][m] that m=0,1,2,3 is for the k_th point in x direction)
		self.coeff = np.zeros((self.n-1, 3, 4)) 
		
		# assigning the values of wrt position and velocity
		self.coeff[:, :, 0] = self.points[0:self.n - 1, :]
		self.coeff[:, :, 1] = self.velo[0:self.n - 1, :] 
		self.coeff[:, :, 2] = 1/self.T*( 3*(self.points[1:self.n, :] - self.points[0:self.n-1, :])/self.T - 2*self.velo[0:self.n-1, :] - self.velo[1:self.n, :])
		self.coeff[:, :, 3] = 1/self.T**2*( 2*(- self.points[1:self.n, :] + self.points[0:self.n-1, :])/self.T + self.velo[0:self.n-1, :] + self.velo[1:self.n, :])

# =================================================================================================
# -- POLYNOMIALS ----------------------------------------------------------------------------------
# =================================================================================================

# here we want to get the coefficient matrix and then build the polynomials accordingly.
# after that we calculate a number of discrete points with this method of interpolation.

class Polynomial:

	def __init__(self, coeff_matrix):
		self.coeff = coeff_matrix

# =================================================================================================
# -- MAIN ----------------------------------------------------------------------------------
# =================================================================================================

print("\n ============================= START ============================= \n")

generator = PositionGenerator(0.3, [0, 0, -0.38], n=100)
cart_position = generator.cart_position()

print("this is our cartesian positions \n\n", cart_position, "\n", type(cart_position), "\n", cart_position.shape, "\n")

print("\n ============================= INVERSE ============================= \n")

n = cart_position.shape[0] # number of points in the cartesian position
THETA = np.zeros((n, 3)) # theta initialization for n points in 3 directions (theta1, theta2, theta3)
for i in range(n):
	inverse = InverseKinematics(cart_position[i])
	inverse.get_J1_positions()
	theta = inverse.get_theta()
	THETA[i, :] = theta

print("this is THETA\n\n", THETA, "\n", type(THETA), "\n", THETA.shape)



	
