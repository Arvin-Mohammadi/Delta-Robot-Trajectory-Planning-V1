# in this file i attempt making a trajectory generator which makes the Delta robot go around a circle 
# the boundary conditions are as followed: 
# 1. continuity on the position 	--> 2n 		conditions
# 2. continuity on the velocity 	--> n - 1 	conditions 
# 3. continuity on the acceleration --> n - 1 	conditions 
# 4. intial and final velocity zero --> 2 		conditions

# =================================================================================================
# -- imports --------------------------------------------------------------------------------------
# =================================================================================================

import time
import numpy as np 
import math 
import matplotlib.pyplot as plt 
from math import sin, cos

# =================================================================================================
# -- FORWARD KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================
# INPUT: actuator angles (unit degrees) 
# OUTPUT: position of EE (unti meters)

def forward_kinematics(Theta, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):

		alpha = np.array(alpha)								# alpha angles
		Theta = np.array(Theta)								# actuator angles in degrees

		# Gearbox ratio = 50, zero offset = 47.2 degree
		Theta = np.array(Theta)/50 - 47.2

		# assigning constants and geometry features 
		theta1 = Theta[0]
		theta2 = Theta[1]
		theta3 = Theta[2]
		e = EE_radius
		f = base_radius
		re = passive_rod
		rf = active_rod

		sqrt3 = 3**0.5
		pi = math.pi
		sin120 = sqrt3/2
		cos120 = -0.5
		tan60 = sqrt3
		sin30 = 0.5
		tan30 = 1/sqrt3

		# forward kinematics: (theta1, theta2, theta3) --> (x0, y0, z0)
		t = (f-e)*tan30/2
		dtr = pi/180
		
		theta1 = theta1*dtr
		theta2 = theta2*dtr
		theta3 = theta3*dtr

		y1 = -(t + rf*math.cos(theta1))
		z1 = -rf*math.sin(theta1)

		y2 = (t + rf*math.cos(theta2))*sin30
		x2 = y2*tan60
		z2 = -rf*math.sin(theta2)

		y3 = (t + rf*math.cos(theta3))*sin30
		x3 = -y3*tan60
		z3 = -rf*math.sin(theta3)

		dnm = (y2-y1)*x3 - (y3-y1)*x2

		w1 = y1**2 + z1**2
		w2 = x2**2 + y2**2 + z2**2
		w3 = x3**2 + y3**2 + z3**2 

		# x = (a1*z + b1)/dnm
		a1 = (z2-z1)*(y3-y1) - (z3-z1)*(y2-y1)
		b1 = -((w2-w1)*(y3-y1) - (w3-w1)*(y2-y1))/2

		# y = (a2*z + b2)/dnm
		a2 = -(z2-z1)*x3 + (z3-z1)*x2
		b2 = ((w2-w1)*x3 - (w3-w1)*x2)/2

		# a*z^2 + b*z + c == 0
		a = a1**2 + a2**2 + dnm**2
		b = 2*(a1*b1 + a2*(b2-y1*dnm) - z1*dnm**2)
		c = (b2 - y1*dnm)*(b2-y1*dnm) + b1**2 + dnm**2*(z1**2 - re**2)

		# discriminant 
		d = b**2 - 4*a*c
		if d < 0:
			return -1

		z0 = -0.5*(b + d**0.5)/a
		x0 = (a1*z0 + b1)/dnm
		y0 = (a2*z0 + b2)/dnm

		return [x0, y0, z0]

# =================================================================================================
# -- INVERSE KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================
# INPUT is EE position (unit meters) and geomety of the Delta robot
# OUTPUT is actuator angles (unit degree)

class InverseKinematics:

	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_triangle_side=0.3464101615, EE_triangle_side=0.2563435195, alpha=[0, 120, 240]):
		
		# initializing the basic geometry and the given data
		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_triangle_side = EE_triangle_side				# the triangle side of the end effoctor e (radius == Triangle_side*(3)**0.5/6)
		self.base_triangle_side = base_triangle_side			# the triangle side of the base or f (radius == Triangle_side*(3)**0.5/6)

	def get_J1_positions(self):
		
		# initializing the J1, F1 position
		self.F1_position = ([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
		J1_position = ([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

		for i in [0, 1, 2]:

			alpha = self.alpha[i]

			# converting the EE position to local position based on alpha
			x = float(self.EE_position_global[0])*cos(alpha*math.pi/180) + float(self.EE_position_global[1])*sin(alpha*math.pi/180)
			y = - float(self.EE_position_global[0])*sin(alpha*math.pi/180) + float(self.EE_position_global[1])*cos(alpha*math.pi/180)
			z = float(self.EE_position_global[2])

			self.EE_position_local = np.array([x, y, z])

			# here positions of the the 3 important points are found (E1, E1_prime, F1) localy 
			tan_30_deg = 0.5773502692
			self.E1_position = self.EE_position_local + np.array([0, -self.EE_triangle_side/2*tan_30_deg, 0]) 
			self.E1_prime_position = np.array([0, self.E1_position[1], self.E1_position[2]])
			self.F1_position[i] = np.array([0, - self.base_triangle_side/2*tan_30_deg, 0])

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
				J1_position = 0
				return

			y = (yF - c1*c2 - c3**0.5)/(c2**2 + 1)
			z = c1 + c2*y

			J1_position[i] = [0, y, z]

		return J1_position

	def get_theta(self, J1_position):

		# finding the angles of actuator joints from the J_1 position and F_1 position
		theta = [0, 0, 0]

		for i in [0, 1, 2]:
			z_J1 = J1_position[i][2]
			y_J1 = J1_position[i][1]
			y_F1 = self.F1_position[i][1]

			theta[i] = (math.atan(-z_J1/(y_F1 - y_J1)))

		# converting to degrees
		theta = np.array(theta)*180/math.pi

		# Gearbox ratio = 50, zero offset = 47.2 degree
		theta = (np.array(theta) + 47.2)*50

		return theta

# =================================================================================================
# -- POSITION GENERATOR ---------------------------------------------------------------------------
# =================================================================================================
# INPUT: ratio (units meters), center cartesian position (units meters), time (units seconds)
# OUTPUT: cartesian points 
class PositionGenerator:

	def __init__(self, ratio, center, t=0.1):
		# t is in seconds
		n = t*1000
		pi = math.pi
		self.n = int(n)
		self.ratio = ratio		# r
		self.center = center	# xc, yc, zc
		self.gamma = np.linspace(0, 2*pi, num=self.n)

	def cartesian_position(self):
		# self.points is the positions of the n points in x, y and z directions (n*3 matrix)
		points = np.array([np.cos(self.gamma)*self.ratio + np.ones((self.n))*self.center[0], np.sin(self.gamma)*self.ratio + np.ones((self.n))*self.center[1], np.ones((self.n))*self.center[2]])
		points = np.transpose(points)
		return points

# =================================================================================================
# -- POLYNOMIAL COEFFICIENT MATRIX ----------------------------------------------------------------
# =================================================================================================
# using the boundary conditions we find a[k][i] for k = 0, ..., n and i = 0, 1, 2, 3

class Coeff:

	def __init__(self, points):

		# initialzing the postion, velocity and acceleration vectors 
		# all of the vectors are n*3 matrices (e.g: for v we have 3 components in x, y and z direction)
		self.points = np.array(points)
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
		velocity_profile = np.zeros(self.points.shape)
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

		velocity_profile[0, :] = initial_velo
		velocity_profile[self.n-1, :] = final_velo
		velocity_profile[1:self.n-1, :] = v

		return velocity_profile

	def coeff_matrix(self, velocity_profile):
		# initializing the coefficient matrix 
		# dim 1 == number of polynomials 				--> k = 0, ..., n-2
		# dim 2 == number of x, y, z directions 		--> 3
		# dim 3 == number of coeff in the polynomial 	--> (e.g: a[k][0][m] that m=0,1,2,3 is for the k_th point in x direction)
		coeff = np.zeros((self.n-1, 3, 4)) 
		
		# assigning the values of wrt position and velocity
		coeff[:, :, 0] = self.points[0:self.n - 1, :]
		coeff[:, :, 1] = velocity_profile[0:self.n - 1, :] 
		coeff[:, :, 2] = 1/self.T*( 3*(self.points[1:self.n, :] - self.points[0:self.n-1, :])/self.T - 2*velocity_profile[0:self.n-1, :] - velocity_profile[1:self.n, :])
		coeff[:, :, 3] = 1/self.T**2*( 2*(- self.points[1:self.n, :] + self.points[0:self.n-1, :])/self.T + velocity_profile[0:self.n-1, :] + velocity_profile[1:self.n, :])
		
		return coeff

# =================================================================================================
# -- POLYNOMIALS ----------------------------------------------------------------------------------
# =================================================================================================
# here we want to get the coefficient matrix and then build the polynomials accordingly.
# after that we calculate a number of discrete points with this method of interpolation.

class Polynomial:

	def __init__(self, coeff_matrix, T):
		
		coeff = coeff_matrix
		self.T = T
		self.n = coeff.shape[0] + 1
		
		# initializing the polynomial matrix
		self.poly_matrix = np.zeros((self.n, 3))
		self.poly_matrix[0, :] = coeff[0, :, 0]
		self.poly_matrix[1:, :] = coeff[:, :, 0] + coeff[:, :, 1]*self.T + coeff[:, :, 2]*self.T**2 + coeff[:, :, 3]*self.T**3
		
		# print("\n this is polynomial matrix\n\n", self.poly_matrix, "\n", type(self.poly_matrix), "\n", self.poly_matrix.shape)

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

generator = PositionGenerator(0.3, [0, 0, -0.38], t=0.1) 	# initiazling position generator 
point_positions = generator.cartesian_position() 	# generating the cartesian coordinates of the points

# calculating IK for the whole profile of movement
theta = np.zeros(point_positions.shape)
for idx, i in enumerate(point_positions):
	inverse = InverseKinematics(i)
	J1_position = inverse.get_J1_positions()
	temp_thata = inverse.get_theta(J1_position)
	theta[idx] = temp_thata

coeff = Coeff(theta) 	# initializing the coeff class (for interpolating theta profile)
velocity_profile = coeff.velocity() 	# calculating velocity profile
coeff_matrix = coeff.coeff_matrix(velocity_profile) 	# calculating coefficient matrix of a_ij
T = coeff.T
t =  coeff.t.transpose()[0]

polynomial = Polynomial(coeff_matrix, T)

plt.grid(True)
plt.plot(t, theta, label=['theta_1', 'theta_2', 'theta_3'])
plt.title("angle-time plot")
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("angle theta (deg)")
plt.savefig("E:\\cubic-spline-1.png")
plt.clf()

plt.grid(True)
plt.plot(t, velocity_profile, label=['theta__dot_1', 'theta__dot_2', 'theta_dot_3'])
plt.title("angular velocity-time plot")
plt.legend()
plt.xlabel("time (ms)")
plt.ylabel("anglular velocity theta_dot (deg/s)")
plt.savefig("E:\\cubic-spline-2.png")
plt.clf()

position = np.zeros(theta.shape)
for idx, i in enumerate(theta):
	position[idx] = forward_kinematics(i)

print(position)

plt.grid(True)
plt.plot(position.transpose()[0], position.transpose()[1])
plt.title("x-y plane")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.savefig("E:\\cubic-spline-3.png")
plt.clf()
