# trajectory planning, point to point movement with 7th order polynomial 
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

# print("time is:")
# print(time.time() - t)

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

		# Gearbox ratio = 50, zero offset = 47.2 degree
		self.theta = ((np.array(self.theta)*180/math.pi) + 47.2)*50

		return self.theta

# =================================================================================================
# -- point to point movement 4-5-6-7 ----------------------------------------------------------------
# =================================================================================================

# in this part we want to move the End-Effector from 1 point to another using a 4-5-6-7 technique 
# in this class the Given Data is: the first and final position of EE, the max value of theta_dot
# and we give the, Desired Data: velocity profile

class PointToPoint4567Movement:
	def __init__(self,EE_position_i, EE_position_f, theta_max=4050):
		self.EE_position_i = EE_position_i
		self.EE_position_f = EE_position_f
		self.theta_max = theta_max*6		# by defaut = 4050 rpm = 4050*6 deg/s
		self.theta_i = np.zeros((3, 1))		# initializing the theta_i 
		self.theta_f = np.zeros((3, 1))		# initializing the theta_f
	
	def inverse_kinematics(self):
		# here we find the Initial and final theta from the Initial and final EE position 

		# Given data of the robot geometry 
		active_rod = 0.2
		passive_rod = 0.46
		base_radius = 0.3464101615 
		EE_radius = 0.2563435195 

		# initializing the InverseKinematics class
		Kinematics = InverseKinematics(self.EE_position_i, active_rod, passive_rod, base_radius, EE_radius)
		Kinematics.get_J1_positions()

		# getting the theta_i from position_i
		theta_i = Kinematics.get_theta()
		
		# getting the theta_f from position_f
		Kinematics.EE_position_global = self.EE_position_f
		Kinematics.get_J1_positions()
		theta_f = Kinematics.get_theta()

		# reshape
		for i in [0, 1, 2]:
			self.theta_i[i] = theta_i[i]
			self.theta_f[i] = theta_f[i]

	def T(self):
		# here we find T(period) from the theta_max
		self.T = 35/16*np.array(self.theta_f - self.theta_i)/self.theta_max
		self.T = max(self.T)

		#print("this is T")
		#print(self.T)

	def theta_dot_t(self):
		# we find the angular velocity profile in this part
		T = math.floor(self.T*1000)
		tau = np.array(range(0, T))/T
		s_tau_d = -140*tau**6 + 420*tau**5 - 420*tau**4 + 140*tau**3
		theta_dot_t = np.array(self.theta_f - self.theta_i)/self.T*s_tau_d

		return theta_dot_t

	def theta_t(self):
		# we find the angular position profile in this part
		T = math.floor(self.T*1000)
		tau = np.array(range(0, T))/T
		s_tau = -20*tau**7 + 70*tau**6 - 84*tau**5 + 35*tau**4
		theta_t = np.array(self.theta_i) + np.array(self.theta_f - self.theta_i)*s_tau

		return theta_t


# =================================================================================================
# -- main -----------------------------------------------------------------------------------------
# =================================================================================================

# initializing the point to point movement class 
movement = PointToPoint4567Movement([0.05, 0.05, -0.31], [0, -0.15, -0.42])
movement.inverse_kinematics()
movement.T()
theta_dot = movement.theta_dot_t()/50 
theta = movement.theta_t()/50 - 47.2

# initializing Position
P = np.zeros(theta.shape)

for i in range(movement.theta_t().shape[1]):
	forward = ForwardKinematics(movement.theta_t()[:, i])
	P[:, i] = forward.get_position()

T = math.floor(movement.T*1000)
tau = np.array(range(0, T))/T


plt.grid(True)
plt.plot(tau, theta.transpose(), label=['theta_1', 'theta_2', 'theta_3'])
plt.title("angle-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("angle theta (deg)")
plt.savefig("E:\\4567-1.png")
plt.clf()

plt.grid(True)
plt.plot(tau, theta_dot.transpose(), label=['theta_dot_1', 'theta_dot_2', 'theta_dot_3'])
plt.title("angular velocity-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("angular velocity theta_dot (deg/s)")
plt.savefig("E:\\4567-2.png")
plt.clf()

plt.grid(True)
plt.plot(tau, P.transpose(), label=['x', 'y', 'z'])
plt.title("position-time plot")
plt.legend()
plt.xlabel("normalized time")
plt.ylabel("position (m)")
plt.savefig("E:\\4567-3.png")
plt.clf()

# converting the data to rpm
# theta_t = movement.theta_t/6
# theta_dot_t = movement.theta_dot_t/6

# writing to file 
# myfile = open("E:\\PointToPoint4567.h", "a")
# list1 = theta_dot_t.tolist()
# line = ["", "", ""]
# for j in [0, 1, 2]:
# 	list1[j] = [round(float(i), 4) for i in list1[j]]
# 	line[j] = ','.join(str(e) for e in list1[j])
# 	myfile.write("float speeds_motor" + str(j+1) + "[] = {" + line[j] + "};\n" )
