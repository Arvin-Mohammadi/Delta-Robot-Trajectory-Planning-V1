# THIS IS NOT USEFUL (and it's not finnished)


# =================================================================================================
# -- imports --------------------------------------------------------------------------------------
# =================================================================================================

import time
import numpy as np 
import math 
import matplotlib.pyplot as plt 
from math import sin, cos
from numpy.linalg import inv 
from numpy import matmul

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
# -- SMOOTHING CUBIC SPLINE -----------------------------------------------------------------------
# =================================================================================================

class SmoothingCubicSpline:

	def __init__(self, points, t_f=0.1):
		# initializing the numpy arrays 
		self.points = np.array(points)
		self.n = self.points.shape[0]-1
		self.t = np.linspace(0, t_f, self.n+1)

	def get_T(self):
		T = self.t[1:self.n+1] - self.t[0:self.n]

		return T

	def get_omega(self, T):
		n = self.n

		# initializing the A matrix 
		A = np.zeros((self.n+1, self.n+1))

		A[0, 0] = 2*T[0]
		A[1, 0] = T[0]

		for i in range(1, self.n):
			A[i-1, i] 	= T[i-1]
			A[i, i]		= 2*(T[i-1] + T[i])
			A[i+1, i]	= T[i]

		A[n-1, n] = T[n-1]
		A[n, n] = 2*T[n-1]


		vi = 0
		vf = 0
		ai = 0
		af = 0
		n = self.n

		# initializing the c matrix 
		C = np.zeros((self.n+1))
		# assigning the values to the C matrix t
		C[0] 	= (self.points[1] 		- self.points[0]	)/T[0] 		- vi
		C[1:n] 	= (self.points[2:n+1] 	- self.points[1:n]	)/T[1:n] 	- (self.points[1:n] - self.points[0:n-1])/T[0:n-1]
		C[n] 	= (self.points[n-1] 	- self.points[n]	)/T[n-1] 	+ vf
		C *= 6

		print(A, A.shape)
		print(C, C.shape)
		return matmul(inv(A), C)

	def get_inv_W(self):
		n = self.n
		W = np.ones((n+1))
		W[0] = 0
		W[n] = 0
		inv_W = np.diag(W)
		# W[0, 0] = 0
		# W[n, n] = 0

		return inv_W

	def get_C(self, T): 
		n = self.n
		C = np.zeros((n+1, n+1))

		C[0, 0] = -1/T[0]
		C[1, 0] = 1/T[0]

		for i in range(1, n):
			C[i-1, i] 	= 1/T[i-1]
			C[i, i]		= -(1/T[i-1] + 1/T[1])
			C[i+1, i]	= 1/T[i]

		C[n-1, n] = 1/T[n-1]
		C[n, n] = -1/T[n-1]

		C *= 6

		return C

	def get_s(self, inv_W, C, omega, moo=0.5):
		n = self.n
		Lambda = (1-moo)/(6*moo) # 0 < moo < 1

		self.points = self.points.reshape((n+1, 1))
		omega = omega.reshape((n+1, 1))
		print(self.points.shape, inv_W.shape, C.shape, omega.shape, sep='\n')
		s = self.points - Lambda*matmul(matmul(inv_W, C), omega)
		return s

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================


# generator = PositionGenerator(0.3, [0, 0, -0.38], t=0.1) 	# initiazling position generator 
# point_positions = generator.cartesian_position() 	# generating the cartesian coordinates of the points



traj = SmoothingCubicSpline([3, -2, -5, 0, 6, 12, 8])
traj.t = np.array([0, 5, 7, 8, 10, 15, 18])
T = traj.get_T()
inv_W = np.diag([0, 1, 1, 1, 1, 1, 0])
C = traj.get_C(T)

omega = traj.get_omega(T)

print(traj.get_s(inv_W, C, omega, moo=0.3))