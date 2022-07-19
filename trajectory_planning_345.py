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

t = time.time()

import numpy as np 
import math 
import matplotlib.pyplot as plt 
print("time is:")
print(time.time() - t)
# =================================================================================================
# -- inverse kinematic class ----------------------------------------------------------------------
# =================================================================================================

# here we have the EE position and the lengths of rods and basic geometry
# our goal is to find theta1, theta2 and theta3 

class InverseKinematics:

	def __init__(self, EE_position, active_rod, passive_rod, base_radius, EE_radius, alpha=[0, 120, 240]):
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
		self.theta = [0, 0, 0]
		for i in [0, 1, 2]:
			z_J1 = self.J1_position[i][2]
			y_J1 = self.J1_position[i][1]
			y_F1 = self.F1_position[i][1]

			self.theta[i] = (math.atan(-z_J1/(y_F1 - y_J1)))
		print(" this is default theta: ")
		print(np.array(self.theta)*180/math.pi)

		# Gearbox ratio = 50, zero offset = 47.2 degree
		self.theta = ((np.array(self.theta)*180/math.pi) + 47.2)*50
		print(" this is theta with GR consideration")
		
		
		return self.theta

# =================================================================================================
# -- point to point movement 3-4-5 ----------------------------------------------------------------
# =================================================================================================

# in this part we want to move the End-Effector from 1 point to another using a 3-4-5 technique 
# in this class the Given Data is: the first and final position of EE, the max value of V
# and we give the, Desired Data: velocity profile

class PointToPoint345Movement:
	def __init__(self,EE_position_i, EE_position_f, theta_max=4050):
		self.EE_position_i = EE_position_i
		self.EE_position_f = EE_position_f
		self.theta_max = theta_max*6		# by defaut = 424.115 rad/s = 4050 rpm
		self.theta_i = np.zeros((3, 1))
		self.theta_f = np.zeros((3, 1))


	def inverse_kinematics(self):

		# here we find the Initial and final theta from the Initial and final EE position 
		active_rod = 0.2
		passive_rod = 0.46
		base_radius = 0.3464101615 
		EE_radius = 0.2563435195 

		Kinematics = InverseKinematics(self.EE_position_i, active_rod, passive_rod, base_radius, EE_radius)
		Kinematics.get_J1_positions()

		theta_i = Kinematics.get_theta()
		
		
		Kinematics.EE_position_global = self.EE_position_f
		Kinematics.get_J1_positions()
		theta_f = Kinematics.get_theta()

		for i in [0, 1, 2]:
			self.theta_i[i] = theta_i[i]
			self.theta_f[i] = theta_f[i]

		print("this is the initial theta")
		print(self.theta_i)
		print(self.theta_i.size)
		print("this is the final theta")
		print(self.theta_f)
		print(self.theta_f.size)

	def T(self):
		# here we find T(period) from the theta_max
		self.T = 15/8*np.array(self.theta_f - self.theta_i)/self.theta_max
		self.T = max(self.T)

		print("this is T")
		print(self.T)


	def theta_dot_t(self):
		T = math.floor(self.T*1000)
		tau = np.array(range(0, T))/T
		s_tau_d = 30*tau**4 - 60*tau**3 + 30*tau**2
		print("this is size of s matrix")
		print(s_tau_d)

		self.theta_dot_t = np.array(self.theta_f - self.theta_i)/self.T*s_tau_d
		print(self.theta_dot_t.shape)
		print("this is theta dot t")
		print(self.theta_dot_t)

	def theta_t(self):
		T = math.floor(self.T*1000)
		tau = np.array(range(0, T))/T
		s_tau = 6*tau**5 - 15*tau**4 + 10*tau**3

		self.theta_t = np.array(self.theta_i) + np.array(self.theta_f - self.theta_i)*s_tau
		print(self.theta_t.shape)
		print("this is theta t")
		print(self.theta_t)

# =================================================================================================
# -- main -----------------------------------------------------------------------------------------
# =================================================================================================


movement = PointToPoint345Movement([0.05, 0.05, -0.31], [0, -0.15, -0.42])
movement.inverse_kinematics()
movement.T()
movement.theta_dot_t()
movement.theta_t()

T = math.floor(movement.T*1000)
tau = np.array(range(0, T))/T
plt.figure()
plt.plot(tau, movement.theta_t.transpose())
plt.figure()
plt.plot(tau, movement.theta_dot_t.transpose())

theta_t = movement.theta_t/6
theta_dot_t = movement.theta_dot_t/6

f = open("E:\\PointToPoint345.h", "a")
f.write("float speeds_motor1[] = { 0")
for i in range(1, T):
	f.write(", " + str(theta_t[0][i])[:7])
f.write("\nfloat speeds_motor2[] = { 0")

for i in range(1, T):
	f.write(", " + str(theta_t[1][i])[:7])
f.write("\nfloat speeds_motor3[] = { 0")

for i in range(1, T):
	f.write(", " + str(theta_t[2][i])[:7])
