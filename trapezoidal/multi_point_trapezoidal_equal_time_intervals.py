
# =================================================================================================
# -- IMPORT ---------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
from math import isclose
import matplotlib.pyplot as plt 
import time 
from scipy.signal import argrelextrema
from math import cos, sin 

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
		# theta = (np.array(theta) + 47.2)*50

		return theta

# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - POINT TO POINT -------------------------------------------------
# =================================================================================================

def ptp_duration(time_profile, k): 				# OUTPUTS: the time difference between two points, acceleration time 
	T = time_profile[-1] - time_profile[0]
	T_a = T*k 
	return T, T_a

def find_acceleration(v_max, T_a):				# OUTPUTS: the acceleration needed for the movement between two points (based on v_max and time)
	return v_max/T_a

# =================================================================================================
# -- POINT TO POINT -------------------------------------------------------------------------------
# =================================================================================================

def	ptp_v_profile(time_interval_profile, v_max, a, T, T_a): 	# velocity profile of the point to point movement 

	velocity_profile = np.zeros(time_interval_profile.shape)

	n_a = 0 



	for idx, t_idx in enumerate(time_interval_profile):
		if t_idx > T_a:
			n_a = idx-1
			break

	velocity_profile[0: n_a] 		= a*time_interval_profile[0: n_a]
	velocity_profile[n_a: -n_a]  	= v_max
	velocity_profile[-n_a:] 		= a*(T - time_interval_profile[-n_a:])
	# velocity_profile[n_a] = 1

	return velocity_profile


# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - MULTI POINT ----------------------------------------------------
# =================================================================================================

def find_time_vector(position):					# OUTPUTS: equally distanced points in time in interval [0, 1] 
	position = np.array(position)
	time_vector = np.linspace(0, 1, position.shape[0])
	return time_vector


def find_v_max(time, position):					# OUTPUTS: the v_max in time intervals [t_i, t_(i+1)]
	time = np.array(time)
	position = np.array(position)
	v_max = np.zeros(time.shape[0])

	for idx, time_idx in enumerate(time[:-1]):
		v_max[idx] = 2*(position[idx+1] - position[idx])/(time[idx+1] - time[idx])

	return v_max 


def find_time_profile(time): 					# OUTPUTS: a time profile based on the time vector [0, T] with the accuracy of 0.01 milisecond
	time = np.array(time)
	time_new = np.arange(0, 1000*time[-1]+1, 1)/1000
	return time_new


def find_idx(time, time_new): 					# OUTPUTS: the index of "time-vector-elements" inside time profile
	time_idx = np.zeros(time.shape)

	counter = 0
	for i, e in enumerate(time_new):
		if e>=time[counter]:
			time_idx[counter] = i
			counter += 1

	return time_idx


def bridge_idx(velocity_profile, v_max): 		# OUTPUTS: the start and finish indeces of the velocity profile points that need to be modified 
	v_max_i 	= v_max[0]
	v_max_ip1 	= v_max[1]
	a = 0 
	b = velocity_profile.shape[0] - 1

	if (v_max_i>0) and (v_max_ip1>0):
		if 	v_max_i <= v_max_ip1:
			for idx, v_idx in enumerate(velocity_profile[:-1]):
				if velocity_profile[idx+1] > v_max_i:
					b = idx
					break
		elif v_max_i > v_max_ip1:
			for idx, v_idx in enumerate(velocity_profile[:-1]):
				if velocity_profile[idx+1] <= v_max_ip1:
					a = idx
					break
	elif (v_max_i<0) and (v_max_ip1<0):
		if 	v_max_i <= v_max_ip1:
			for idx, v_idx in enumerate(velocity_profile[:-1]):
				if velocity_profile[idx+1] > v_max_ip1:
					a = idx
					break
		elif v_max_i > v_max_ip1:
			for idx, v_idx in enumerate(velocity_profile[:-1]):
				if velocity_profile[idx+1] < v_max_i:
					b = idx
					break
	return (a, b)

# =================================================================================================
# -- MULTI POINT ----------------------------------------------------------------------------------
# =================================================================================================

def find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile): 	# OUTPUTS: the velocity profile and position profile
	time_index = find_idx(time_vector, (time_profile))
	velocity_profile = np.zeros(time_profile.shape)
	position_profile = np.zeros(time_profile.shape)


	k = 0.5 												# k = T_a/T
	
	for idx, time_idx in enumerate(time_vector[:-1]):
		idx_i 	= int(time_index[idx])
		idx_ip1 = int(time_index[idx+1])

		# defining time intervals
		time_interval_profile = np.copy(time_profile[idx_i:idx_ip1+1])
		time_interval_profile -= time_interval_profile[0] 

		v_max = v_max_vector[idx]							# v_max with in time interval 

		T, T_a 	= ptp_duration(time_vector[idx:idx+2], k) 	# duration of point to point movement from start to finish 
		a 		= find_acceleration(v_max, T_a) 			# returns the measure of acceleration within time interval
		
		velocity_profile[idx_i:idx_ip1+1] = ptp_v_profile(time_interval_profile, v_max, a, T, T_a)

		velocity_profile[idx_i] = 0
		velocity_profile[idx_ip1] = 0


	time_step = 1e-03
	current_pos = position_vector[0] 						# calculating the integral of the velocity profile = position
	for idx, v in enumerate(velocity_profile[1:]):
		current_pos += v*time_step
		position_profile[idx+1] = current_pos

	return (velocity_profile, position_profile)


def modify_velocity_profile(time_profile, velocity_profile, v_max_vector, p0): 		# OUTPUTS: the new modified velocity and position profile (with the bridge applied)

	position_profile = np.zeros(time_profile.shape)

	v_max_idx_vector = np.zeros(v_max_vector.shape) # index of "v_max_vector-elements" in v profile
	
	counter = 0
	for idx, v_idx in enumerate(velocity_profile):
		if abs(v_idx - v_max_vector[counter]) < 10**(-10):
			v_max_idx_vector[counter] = int(idx)
			counter += 1


	for idx, v_max_idx in enumerate(v_max_vector[:-1]):
		a, b = bridge_idx(velocity_profile[int(v_max_idx_vector[idx]):int(v_max_idx_vector[idx+1])+3], v_max_vector[idx:idx+2])
		a += int(v_max_idx_vector[idx])
		b += int(v_max_idx_vector[idx])

		if (v_max_vector[idx]*v_max_vector[idx+1]>0):
			velocity_profile[a:b+1] = np.sign(v_max_vector[idx])*min(abs(v_max_vector[idx]), abs(v_max_vector[idx+1]))
		else:
			pass 

	time_step = 1e-03
	current_pos = p0
	for idx, v in enumerate(velocity_profile[1:]):
		current_pos += v*time_step
		position_profile[idx+1] = current_pos

	return (velocity_profile, position_profile)


def finalized_velocity_profile(time_profile, velocity_profile, v_max_vector, p0): 		# OUTPUTS: finalized velocity and position profile (with the bridge shortened)

	position_profile = np.zeros(time_profile.shape)
	velocity_profile_old = np.copy(velocity_profile)
	check = -1
	a = -1 
	b = -1 

	for i, v_i in enumerate(velocity_profile[:-1]):
		
		if i <= check:
			continue

		if (velocity_profile[i] == velocity_profile[i+1]) and (velocity_profile[i] != 0):
			a = i
			for j, v_j in enumerate(velocity_profile[i+1:]):
				if v_j != velocity_profile[i]:
					check = i + j + 1
					b = i + j
					break

			if abs(int(a) - int(b)) <= 3:
				continue 
			else: 
				temp_len = velocity_profile[b: ].shape[0]
				starting_idx = int(np.floor((b+a+1)/2))

				velocity_profile[starting_idx: starting_idx+temp_len] = velocity_profile[b:]
				velocity_profile[starting_idx+temp_len-2:] = 0

	time_step = 1e-03
	current_pos = p0
	for idx, v in enumerate(velocity_profile[1:]):
		current_pos += v*time_step
		position_profile[idx+1] = current_pos

	return (velocity_profile, position_profile)

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

def main(position_vector):

	# the vectors 
	time_vector = find_time_vector(position_vector)			# define the time vector (equally distanced times)
	v_max_vector = find_v_max(time_vector, position_vector)	# calculate v_max in each time interval

	# porfiles
	time_profile = find_time_profile(time_vector)																												# returns (new time profile (miliseconds), index of the time vector elements in the time profile)
	velocity_profile, position_profile = find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile)										# vector profile (trapezoidal)

	velocity_profile_new, position_profile_new = modify_velocity_profile(time_profile, np.copy(velocity_profile), v_max_vector, position_vector[0]) 			# the new velocity and position profile with the bridge applied

	velocity_profile_final, position_profile_final = finalized_velocity_profile(time_profile, np.copy(velocity_profile_new), v_max_vector, position_vector[0]) 	# the final velocity and position profile with the bridge shortened

	plt.plot(time_profile, velocity_profile)
	plt.plot(time_profile, velocity_profile_new)
	plt.plot(time_profile, velocity_profile_final)
	plt.show()
	plt.plot(time_profile, position_profile)
	plt.plot(time_profile, position_profile_new)
	plt.plot(time_profile, position_profile_final)
	plt.plot(time_vector, position_vector, 'ro')
	plt.show()



# =================================================================================================
# -- TEST 1 ---------------------------------------------------------------------------------------
# =================================================================================================

# main([0, 0.05])
main([0, 0.02, 0.05, 0.0])
main([0, 0.02, 0.03, 0.05, -0.03, 0.0])
# main([0, 0.02, 0.05, 0.04, -0.02, -0.08, 0.0])
# main([0, 0.02, 0.05, 0.04, -0.02, -0.06, -0.08, -0.05, 0.0])

# =================================================================================================
# -- TEST 2 ---------------------------------------------------------------------------------------
# =================================================================================================

# theta = np.linspace(0, 2*np.pi, 100)
# x = np.cos(theta)*0.30
# y = np.sin(theta)*0.30
# z = -0.35

# angle1 = []
# for i in range(100):
# 	inverse = InverseKinematics([x[i], y[i], z])

# 	J1 = inverse.get_J1_positions()

# 	temp_theta = inverse.get_theta(J1)
# 	angle1.append(temp_theta[0])

# main(x)
