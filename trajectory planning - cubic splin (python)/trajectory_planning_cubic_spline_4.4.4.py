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
# -- FORWARD KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================

# here we have the angles (theta1, theta2, theta3) and the geometry of the robot 
# our goal is to find EE position based on the given data


class ForwardKinematics:

	def __init__(self, Theta, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):
		# initializing the basic geometry and the given data

		self.alpha = np.array(alpha)							# alpha angles
		self.Theta = np.array(Theta)/50 - 47.2					# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_radius = EE_radius								# the radius of the end effoctor e 
		self.base_radius = base_radius							# the radius of the base or f
	
	def get_position(self):
		# assigning constants and geometry features 
		theta1 = self.Theta[0]
		theta2 = self.Theta[1]
		theta3 = self.Theta[2]
		e = self.EE_radius
		f = self.base_radius
		re = self.passive_rod
		rf = self.active_rod

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
		# self.theta = ((np.array(self.theta)*180/math.pi) + 47.2)*50

		return self.theta

# =================================================================================================
# -- POSITION GENERATOR ---------------------------------------------------------------------------
# =================================================================================================

# in this part we try to generate a circle and get n sample points from that circle. 
# this results in having n-1 cartesian points (matrix.shape=(n-1, 3))

class PositionGenerator:

	def __init__(self, ratio, center, t=6):
		# t is in seconds
		n = 500
		pi = math.pi
		self.n = int(n)
		self.ratio = ratio		# r
		self.center = center	# xc, yc, zc
		self.gamma = np.linspace(0, 2*pi, num=self.n+1)

	def cart_position(self):
		# self.points is the positions of the n points in x, y and z directions (n*3 matrix)
		self.points = np.array([np.cos(self.gamma)*self.ratio + np.ones((self.n+1))*self.center[0], np.sin(self.gamma)*self.ratio + np.ones((self.n+1))*self.center[1], np.ones((self.n+1))*self.center[2]])
		self.points = np.transpose(self.points)
		return self.points

# =================================================================================================
# -- POLYNOMIAL COEFFICIENT MATRIX ----------------------------------------------------------------
# =================================================================================================

# Give data: n-1 points

# first step: making the n-1 points into n+1 points (adding q_1 and q_(n-1))
# second step: 

class Coeff:

	def __init__(self, points, final_time, intial_velo=[0, 0, 0], final_velo=[0, 0, 0], initial_acc=[0, 0, 0], final_acc=[0, 0, 0]):

		# initialzing the postion, velocity and acceleration vectors 
		# all of the vectors are n*3 matrices (e.g: for v we have 3 components in x, y and z direction)
		self.points = np.array(points)
		self.t = np.zeros(self.points.shape)
		self.n = self.points.shape[0] - 1
		self.vi = np.array(intial_velo)
		self.vf = np.array(final_velo)
		self.ai = np.array(initial_acc)
		self.af = np.array(final_acc)

		# time value assignment, shape = n*3 it is measured in milisecond
		for i in [0, 1, 2]:
			self.t[:, i] = np.linspace(0, final_time, num=self.n+1)

		# initializing q-bar matrix
		self.q_bar = np.zeros((self.points.shape[0]+2, self.points.shape[1]))
		# assigning values to q-bar matrix without assigning q_1 and q_n_minut_1
		self.q_bar[0, :] = self.points[0, :]
		self.q_bar[2:self.n-1, :] = self.points[1:self.n-2]
		self.q_bar[self.n, :] = self.points[self.n-2, :]
		
	def get_t_bar(self):
		# initializing t-bar matrix  
		self.t_bar = np.zeros((self.t.shape[0]+2, self.t.shape[1]))
		# assigning values to t-bar matrix
		self.t_bar[0, :] = self.t[0, :]
		self.t_bar[2:self.n-1, :] = self.t[1:self.n-2]
		self.t_bar[self.n, :] = self.t[self.n-2, :]
		t_1 = (self.t[0, :] + self.t[1, :])/2
		t_n_minus_1 = (self.t[self.n-3, :] + self.t[self.n-2, :])/2
		self.t_bar[[1, self.n-1], :] = [t_1, t_n_minus_1]
	
	def get_T(self):
		# intializing the T matrix (time intervals)
		self.T = np.zeros((self.n, 3))
		# assigning values to T matrix 
		self.T = self.t_bar[1:self.n+1, :] - self.t_bar[0:self.n, :] # changed
	
	def get_omega(self):
		
		vi = self.vi
		vf = self.vf
		ai = self.ai
		af = self.af
		n = self.n

		# initializing the c matrix 
		c = np.zeros((self.n-1, 3))
		# assigning the values to the c matrix t
		c[0, :] = (self.q_bar[2] - self.q_bar[0])/self.T[1] - vi*(1 + self.T[0]/self.T[1]) - ai*(0.5 + self.T[0]/(3*self.T[1]))*self.T[0]
		c[1, :] = (self.q_bar[3] - self.q_bar[2])/self.T[2] - (self.q_bar[2] - self.q_bar[0])/self.T[1] + vi*self.T[0]/self.T[1] + ai*self.T[0]**2/(3*self.T[1])
		c[2:n-3, :] = (self.q_bar[4:n-1] - self.q_bar[3:n-2])/self.T[3:n-2] - (self.q_bar[3:n-2] - self.q_bar[2:n-3])/self.T[2:n-3]
		c[n-3, :] = (self.q_bar[n] - self.q_bar[n-2])/self.T[n-2] - (self.q_bar[n-2] - self.q_bar[n-3])/self.T[n-3] - vf*self.T[n-1]/self.T[n-2] + af*self.T[n-1]**2/(3*self.T[n-2])
		c[n-2, :] =(self.q_bar[n-2] - self.q_bar[n])/self.T[n-2] + vf*(1 + self.T[n-1]/self.T[n-2]) - af*(0.5 + self.T[n-1]/(3*self.T[n-2])*self.T[n-1])
		c = c*6

		# initializing the A matrix 
		A = np.zeros((3, self.n-1,self.n-1))
		# assigning the values to the A matrix 
		for i in range(1, A.shape[1]-1): # changed
			A[:, i-1, i] = self.T[i]
			A[:, i, i] = 2*(self.T[i] + self.T[i+1])
			A[:, i+1, i] = self.T[i+1]
		
		
		A[:, 0, 0] = 2*self.T[1] + self.T[0]*(3 + self.T[0]/self.T[1])
		A[:, 1, 0] = self.T[1] - self.T[0]**2/self.T[1]
		A[:, n-3, n-2] = self.T[n-2] - self.T[n-1]**2/self.T[n-2]
		A[:, n-2, n-2] = 2*self.T[n-2] + self.T[n-1]*(3 + self.T[n-1]/self.T[n-2])

		
		# initializing omega matrix 
		self.omega = np.zeros((self.n-1, 3))
		for i in [0, 1, 2]:
			M = np.linalg.inv(A[i, :, :])
			N = c[:, i]
			self.omega[:, i] = np.matmul(M, N)

	def get_q_bar(self):
		
		vi = self.vi
		vf = self.vf
		ai = self.ai
		af = self.af
		n = self.n
		# updating the q bar index 1 an (n-1) with the help of omega matrix 
		self.q_bar[1, :] = self.q_bar[0] + self.T[0]*vi + self.T[0]**2/3*ai + self.T[0]**2/6*self.omega[0]
		self.q_bar[n-1, :] = self.q_bar[n] - self.T[n-1]*vf + self.T[n-1]**2/3*af + self.T[n-1]**2/6*self.omega[n-2]

	def get_coeff_matrix(self):
		# initializing the coefficient matrix 
		# dim 1 == number of polynomials 				--> k = 0, ..., n-1
		# dim 2 == number of x, y, z directions 		--> 3
		# dim 3 == number of coeff in the polynomial 	--> (e.g: a[k][0][m] that m=0,1,2,3 is for the k_th point in x direction)
		self.coeff = np.zeros((self.n, 3, 4)) 
		
		# omega matrix has n-1 rows but it needs n+1 so we re-assign the omega to a n+1*3 matrix putting the index 0 and n to zero
		
		temp_omega = self.omega
		self.omega = np.zeros((self.n+1, 3))
		self.omega[0, :] = 0
		self.omega[self.n, :] = 0
		self.omega[1:self.n, :] = temp_omega

		# assigning the values of wrt position and velocity
		self.coeff[:, :, 0] = self.q_bar[0:self.n, :]
		self.coeff[:, :, 1] = (self.q_bar[1:self.n+1, :] - self.q_bar[0:self.n, :])/self.T[0:self.n, :] - self.T[0:self.n, :]/6*(self.omega[1:self.n+1, :] + 2*self.omega[0:self.n, :])
		self.coeff[:, :, 2] = self.omega[0:self.n, :]/2
		self.coeff[:, :, 3] = (self.omega[1:self.n+1, :] - self.omega[0:self.n, :])/(6*self.T[0:self.n, :])
		
		return self.coeff

# =================================================================================================
# -- POLYNOMIALS ----------------------------------------------------------------------------------
# =================================================================================================

# here we want to get the coefficient matrix and then build the polynomials accordingly.
# after that we calculate a number of discrete points with this method of interpolation.

class Polynomial:

	def __init__(self, coeff_matrix, T):
		
		self.coeff = coeff_matrix
		self.T = T
		self.n = self.coeff.shape[0]
		
		# initializing the polynomial matrix
		self.poly_matrix = np.zeros((self.n+1, 3))

		self.poly_matrix[0, :] = self.coeff[0, :, 0]
		self.poly_matrix[1:, :] = self.coeff[:, :, 0] + self.coeff[:, :, 1]*self.T + self.coeff[:, :, 2]*self.T**2 + self.coeff[:, :, 3]*self.T**3
		
		print("\n this is polynomial theta matrix\n\n", self.poly_matrix, "\n", type(self.poly_matrix), "\n", self.poly_matrix.shape)
	
	def get_velo(self):
		self.velocity = np.zeros((self.n+1, 3))
		self.velocity[0, :] = self.coeff[0, :, 1]
		self.velocity[1:, :] = self.coeff[:, :, 1] + 2*self.coeff[:, :, 2]*self.T + 3*self.coeff[:, :, 3]*self.T**2
		
		
		print("\n this is polynomial angular velocity matrix\n\n", self.velocity, "\n", type(self.velocity), "\n", self.velocity.shape)
		print("maximum point is: ", np.argmax(self.velocity))

		return self.velocity 

	def get_acc(self):
		self.acceleration = np.zeros((self.n+1, 3))
		self.velocity[0, :] = 2*self.coeff[0, :, 2]
		self.acceleration[1:, :] = 2*self.coeff[:, :, 2] + 6*self.coeff[:, :, 3]*self.T

		print("\n this is polynomial angular acceleration matrix\n\n", self.acceleration, "\n", type(self.acceleration), "\n", self.acceleration.shape)

		return self.acceleration


# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

# print("\n ============================= START ============================= \n")

generator = PositionGenerator(0.3, [0, 0, -0.38], t=6)
cart_position = generator.cart_position()

#print("this is our cartesian positions \n\n", cart_position, "\n", type(cart_position), "\n", cart_position.shape, "\n")

#print("\n ============================= INVERSE ============================= \n")

n = cart_position.shape[0] # number of points in the cartesian position
THETA = np.zeros((n, 3)) # theta initialization for n points in 3 directions (theta1, theta2, theta3)
for i in range(n):
	inverse = InverseKinematics(cart_position[i])
	inverse.get_J1_positions()
	theta = inverse.get_theta()
	THETA[i, :] = theta

#print("this is THETA\n\n", THETA, "\n", type(THETA), "\n", THETA.shape)

#print("\n ============================= POLY COEFF ============================= \n")

coeff = Coeff(THETA, 6)
coeff.get_t_bar()
coeff.get_T()
coeff.get_omega()
coeff.get_q_bar()
coeff_matrix = coeff.get_coeff_matrix()
t = coeff.t
T = coeff.T

# print("\n this is polynomial coefficients matrix\n\n", coeff_matrix, "\n", type(coeff_matrix), "\n", coeff_matrix.shape)

print("\n ============================= POLY DESCRETE POINTS ============================= \n")

polynomial = Polynomial(coeff_matrix, T)
VELOCITY = polynomial.get_velo()
ACCELERATION = polynomial.get_acc()

print("time is:")
print(time.time() - t)

# print("\n ============================= FORWARD KINEMATICS ============================= \n")

POSITION = np.zeros((n, 3)) # theta initialization for n points in 3 directions (theta1, theta2, theta3)
for i in range(n):
	forward = ForwardKinematics(THETA[i])
	position = forward.get_position()
	POSITION[i, :] = position


print("this is t\n", t, "\n", t.shape)
# print("\n this is cartesian position using forward kinetmatics\n\n", POSITION, "\n", type(POSITION), "\n", POSITION.shape)

plt.title("position in x-y plane ", fontsize='16')
plt.plot(POSITION[:, 0], POSITION[:, 1])
plt.xlabel("X", fontsize='13')
plt.ylabel("Y", fontsize='13')
plt.legend(('POSITION'), loc='best')
plt.grid()
plt.show()

plt.title("velocity in x, y and z directions", fontsize='16')
print("hey i'm HERE!!!", VELOCITY[:, 0].shape)
plt.plot(t[:, 0], VELOCITY[:, 0]/6, label="1")
plt.plot(t[:, 1], VELOCITY[:, 1]/6, label="2")
plt.plot(t[:, 2], VELOCITY[:, 2]/6, label="3")
plt.legend()
plt.xlabel("t")
plt.ylabel("velocity")
plt.grid()
plt.show()

plt.title("acceleration in x, y and z directions", fontsize='16')
plt.plot(t, ACCELERATION[:, 0])
plt.plot(t, ACCELERATION[:, 1])
plt.plot(t, ACCELERATION[:, 2])
plt.xlabel("t")
plt.ylabel("acceleration")
plt.grid()
plt.show()

# # =================================================================================================
# # -- TEST -----------------------------------------------------------------------------------------
# # =================================================================================================

inverse = InverseKinematics([0.0, 0.0, -0.38])
inverse.get_J1_positions()
print(inverse.get_theta()*180/math.pi)

print("time is:")
print(time.time() - t)
