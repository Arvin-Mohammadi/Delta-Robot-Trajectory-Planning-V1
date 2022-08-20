# attempt #1 
# trajectory planning, point to point movement with 5rd order polynomial 
# 2 fundamental questions: 
# 1. given a specific location in the real world, what values should my robot's joint be ...
# ... set to in order to get the EE there? (inverse kinematics)
# 2. given the setting of my joints, where is my EE in real world coordinates? (forward kinematics)

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
# -- point to point movement 3-4-5 ----------------------------------------------------------------
# =================================================================================================

# in this part we want to move the End-Effector from 1 point to another using a 3-4-5 technique 
# in this class the Given Data is: the first and final position of EE, the max value of V
# and we give the, Desired Data: velocity profile

class PointToPoint345Movement:
	def __init__(self,EE_position_i, EE_position_f, theta_dot_max=4050):
		self.EE_position_i = EE_position_i
		self.EE_position_f = EE_position_f
		self.theta_dot_max = theta_dot_max*6		# by defaut = 424.115 rad/s = 4050 rpm
		theta_i = np.zeros((3, 1))
		theta_f = np.zeros((3, 1))


	def inverse_kinematics(self):

		Kinematics = InverseKinematics(self.EE_position_i)
		J1_position = Kinematics.get_J1_positions()
		theta_i = Kinematics.get_theta(J1_position)
		
		
		Kinematics = InverseKinematics(self.EE_position_f)
		J1_position = Kinematics.get_J1_positions()
		theta_f = Kinematics.get_theta(J1_position)

		theta_i = theta_i.reshape((3, 1)) 
		theta_f = theta_f.reshape((3, 1))

		return (theta_i, theta_f)


	def T(self, theta_i, theta_f):
		# here we find T(period) from the theta_dot_max
		T = 15/8*np.array(theta_f - theta_i)/self.theta_dot_max
		T = max(T)
		
		return T

	def theta_dot_t(self, theta_i, theta_f, T):
		# we find the angular velocity profile in this part
		T = math.floor(T*1000)
		tau = np.array(range(0, T))/T
		s_tau_d = 30*tau**4 - 60*tau**3 + 30*tau**2
		# print("this is size of s matrix")
		# print(s_tau_d)

		theta_dot_t = np.array(theta_f - theta_i)/T*s_tau_d

		return theta_dot_t

	def theta_t(self, theta_i, theta_f, T):
		# we find the angular position profile in this part
		T = math.floor(T*1000)
		tau = np.array(range(0, T))/T
		s_tau = 6*tau**5 - 15*tau**4 + 10*tau**3

		theta_t = np.array(theta_i) + np.array(theta_f - theta_i)*s_tau

		return theta_t

# =================================================================================================
# -- main -----------------------------------------------------------------------------------------
# =================================================================================================

movement = PointToPoint345Movement([0.05, 0.05, -0.31], [0, -0.15, -0.42])	# initizalizing the movement 
theta_i, theta_f = movement.inverse_kinematics()	# calculating the ineverse kinemtatics for initial and final points 

T = movement.T(theta_i, theta_f)	# calculating T 
theta = movement.theta_t(theta_i, theta_f, T)	# calculaing theta profile 
theta_dot = movement.theta_dot_t(theta_i, theta_f, T) # calculating theta_dot profile

T = math.floor(T*1000)
tau = np.array(range(0, T))/T 	# normalized time 

position = np.zeros(theta.shape)
for idx, i in enumerate(theta.transpose()):
	print(idx)
	position[:, idx] = forward_kinematics(theta[:, idx])

plt.grid(True)
plt.plot(tau, theta.transpose(), label=['theta_1', 'theta_2', 'theta_3'])
plt.title("angle-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("angle theta (deg)")
plt.savefig("E:\\345-1.png")
plt.clf()

plt.grid(True)
plt.plot(tau, theta_dot.transpose(), label=['theta__dot_1', 'theta__dot_2', 'theta__dot_3'])
plt.title("angular velocity-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("angular velocity theta_dot (deg/s)")
plt.savefig("E:\\345-2.png")
plt.clf()

plt.grid(True)
plt.plot(tau, position.transpose(), label=['x', 'y', 'z'])
plt.title("position-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("position (m/s)")
plt.savefig("E:\\345-3.png")
plt.clf()

# == writing the output file =======================================================================
# f = open("E:\\PointToPoint345.h", "a")
# f.write("float speeds_motor1[] = { 0")
# for i in range(1, T):
# 	f.write(", " + str(theta_t[0][i])[:7])
# f.write("\nfloat speeds_motor2[] = { 0")

# for i in range(1, T):
# 	f.write(", " + str(theta_t[1][i])[:7])
# f.write("\nfloat speeds_motor3[] = { 0")

# for i in range(1, T):
# 	f.write(", " + str(theta_t[2][i])[:7])

