# in this file i attempt making a trajectory generator which makes the Delta robot go around a circle 
# the boundary conditions are as followed: (we have n-1 points in the trajectory and then we add another 2 ponits)
# 1. continuity on the position 		 --> 2n-2 	conditions
# 2. continuity on the velocity 		 --> n-2 	conditions 
# 3. continuity on the acceleration 	 --> n-2 	conditions 
# 4. initial and final velocity zero 	 --> 2 		conditions
# 5. initial and final acceleration zero --> 2 		conditions

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

import time
t = time.time()

import numpy as np 
import math 
import matplotlib.pyplot as plt 

print("time is:")
print(time.time() - t)


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
				print('\n non existing point \n')
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
# this results in having n-1 cartesian points (matrix.shape=(n-1, 3))

class PositionGenerator:

	def __init__(self, ratio, center, t=0.1):
		# t is in seconds
		n = t*1000
		pi = math.pi
		self.n = int(n)
		self.ratio = ratio		# r
		self.center = center	# xc, yc, zc
		self.gamma = np.linspace(0, 2*pi, num=self.n-1)

	def cart_position(self):
		# self.points is the positions of the n points in x, y and z directions (n*3 matrix)
		self.points = np.array([np.cos(self.gamma)*self.ratio + np.ones((self.n-1))*self.center[0], np.sin(self.gamma)*self.ratio + np.ones((self.n-1))*self.center[1], np.ones((self.n-1))*self.center[2]])
		self.points = np.transpose(self.points)
		return self.points

# =================================================================================================
# -- POLYNOMIAL COEFFICIENT MATRIX ----------------------------------------------------------------
# =================================================================================================

# Give data: n-1 points

# first step: making the n-1 points into n+1 points (adding q_1 and q_(n-1))
# second step: 

class Coeff:

	def __init__(self, points, intial_velo=[0, 0, 0], final_velo=[0, 0, 0], initial_acc=[0, 0, 0], final_acc=[0, 0, 0]):

		# initialzing the postion, velocity and acceleration vectors 
		# all of the vectors are n*3 matrices (e.g: for v we have 3 components in x, y and z direction)
		self.points = np.array(points)
		self.t = np.zeros(self.points.shape)
		self.n = self.points.shape[0] + 1
		self.vi = np.array(intial_velo)
		self.vf = np.array(final_velo)
		self.ai = np.array(initial_acc)
		self.af = np.array(final_acc)

		# time value assignment, shape = n*3 it is measured in milisecond
		for i in [0, 1, 2]:
			self.t[:, i] = np.linspace(0, self.n, num=self.n-1)
		
		# # initializing t-bar matrix  
		# self.t_bar = np.zeros((self.t.shape[0]+2, self.t.shape[1]))
		# # assigning values to t-bar matrix
		# self.t_bar[0, :] = self.t[0, :]
		# self.t_bar[2:self.n-1, :] = self.t[1:self.n-2]
		# self.t_bar[self.n, :] = self.t[self.n-2, :]
		# t_1 = (self.t[0, :] + self.t[1, :])/2
		# t_n_minus_1 = (self.t[self.n-3, :] + self.t[self.n-2, :])/2
		# self.t_bar[[1, self.n-1], :] = [t_1, t_n_minus_1]

		# intializing the T matrix (time intervals)
		self.T = np.zeros((self.n, 3))
		# assigning values to T matrix 
		self.T = self.t[1:self.n-1, :] - self.t[0:self.n-2, :] # changed

		# # initializing q-bar matrix
		# self.q_bar = np.zeros((self.points.shape[0]+2, self.points.shape[1]))
		# # assigning values to q-bar matrix with out assigning q_1 and q_n_minut_1
		# self.q_bar[0, :] = self.points[0, :]
		# self.q_bar[2:self.n-1, :] = self.points[1:self.n-2]
		# self.q_bar[self.n, :] = self.points[self.n-2, :]

	def get_omega(self):
		# initializing the c matrix 
		c = np.zeros((self.n-1, 3))
		# assigning the values to the c matrix 
		c[0, :] = (self.points[2, :] - self.points[0, :])/self.T[1, :] - self.vi*(1 + self.T[0, :]/self.T[1, :]) - self.ai*(0.5 + self.T[0, :]/(3*self.T[1, :]))*self.T[0, :]
		c[1, :] = (self.points[3, :] - self.points[2, :])/self.T[2, :] - (self.points[2, :] - self.points[0, :])/self.T[1, :] + self.vi*self.T[0, :]/self.T[1, :] + self.ai*self.T[0, :]**2/(3*self.T[1, :])
		c[2:self.n-3, :] = (self.points[4:self.n-1, :] - self.points[3:self.n-2, :])/self.T[3:self.n-2, :] - (self.points[3:self.n-2, :] - self.points[2:self.n-3, :])/self.T[2:self.n-3, :]
		c[self.n-3, :] = (self.points[self.n-2, :] - self.points[self.n-4, :])/self.T[self.n-4, :] - (self.points[self.n-4, :] - self.points[self.n-5, :])/self.T[self.n-5, :] - self.vf*self.T[self.n-3, :]/self.T[self.n-4, :] + self.af*self.T[self.n-3, :]**2/(3*self.T[self.n-4, :]) # changed
		c[self.n-2, :] = ( - self.points[self.n-2, :] + self.points[self.n-4, :])/self.T[self.n-4, :] + self.vf*(1 + self.T[self.n-3 ,:]/self.T[self.n-4, :]) - self.af*(0.5 + self.T[self.n-3, :]/(3*self.T[self.n-4]))*self.T[self.n-3, :] # changed 
		c = c*6

		# initializing the A matrix 
		A = np.zeros((3, self.n-1,self.n-1))
		# assigning the values to the A matrix 
		for i in range(1, A.shape[1]-3): # changed
			A[:, i-1, i] = self.T[i]
			A[:, i, i] = 2*(self.T[i] + self.T[i+1])
			A[:, i+1, i] = self.T[i+1]
		A[:, 0, 0] = 2*self.T[1] + self.T[0]*(3 + self.T[0]/self.T[1])
		A[:, 1, 0] = self.T[1] - self.T[0]**2/self.T[2]
		A[:, self.n-4, self.n-4] = 2*self.T[self.n-4] + self.T[self.n-3]*(3 + self.T[self.n-3]/self.T[self.n-4]) # changed
		A[:, self.n-5, self.n-4] = self.T[self.n-4] - self.T[self.n-3]**2/self.T[self.n-4]	# changed 

		print("\nthis is A\n\n")
		print(A)
		print(A.shape)
		print("\n this is c matrix\n\n")
		print(c)
		print(c.shape)
		# initializing omega matrix 
		self.omega = np.zeros((self.n-1, 3))
		for i in [0, 1, 2]:
			M = np.linalg.inv(A[i, :, :])
			N = c[:, i]
			self.omega[:, i] = np.matmul(M, N)
		
		
		print("\n this is omega \n\n")
		print(self.omega)
		print(self.omega.shape)
		return self.omega
	
# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

print("\n ============================= START ============================= \n")

generator = PositionGenerator(0.3, [0, 0, -0.38], t=0.1)
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

print("\n ============================= POLY COEFF ============================= \n")

#coeff = Coeff(THETA)
#coeff.get_omega()
#print("\n this is polynomial coefficients matrix\n\n", coeff_matrix, "\n", type(coeff_matrix), "\n", coeff_matrix.shape)
	
# =================================================================================================
# -- TEST -----------------------------------------------------------------------------------------
# =================================================================================================

print("\n ============================= TEST ============================= \n")

points = [[3, 0, 0], [-2, 0, 0], [-5, 0, 0], [0, 0, 0], [6, 0, 0], [12, 0, 0], [8, 0, 0]]

coeff = Coeff(points)
coeff.t = [[0, 0, 0], [2.5, 2.5, 2.5], [5, 5, 5], [7, 7, 7], [8, 8, 8], [10, 10, 10], [15, 15, 15], [16.5, 16.5, 16.5], [18, 18, 18]]
coeff.t = np.array(coeff.t)
coeff.T = [[2.5, 2.5, 2.5], [2.5, 2.5, 2.5], [2, 2, 2], [1, 1, 1], [2, 2, 2], [5, 5, 5], [1.5, 1.5, 1.5], [1.5, 1.5, 1.5]]
coeff.T = np.array(coeff.T)

coeff.get_omega()