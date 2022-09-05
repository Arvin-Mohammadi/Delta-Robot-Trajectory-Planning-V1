
# =================================================================================================
# -- IMPORT ---------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
from math import isclose
import matplotlib.pyplot as plt 
import time 
from scipy.signal import argrelextrema

# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - POINT TO POINT -------------------------------------------------
# =================================================================================================

def ptp_duration(time_profile, k):
	T = time_profile[-1] - time_profile[0]
	T_a = T*k 
	return T, T_a

def find_acceleration(v_max, T_a):
	return v_max/T_a

# =================================================================================================
# -- POINT TO POINT -------------------------------------------------------------------------------
# =================================================================================================

def	ptp_v_profile(time_interval_profile, v_max, a, T, T_a):

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

def find_time_vector(position):	# OUTPUTS: equally distanced points in time in interval [0, 1] 
	position = np.array(position)
	time_vector = np.linspace(0, 1, position.shape[0])
	return time_vector


def find_v_max(time, position):	# OUTPUTS: the v_max in time intervals [t_i, t_(i+1)]
	time = np.array(time)
	position = np.array(position)
	v_max = np.zeros(time.shape[0])

	for idx, time_idx in enumerate(time[:-1]):
		v_max[idx] = 2*(position[idx+1] - position[idx])/(time[idx+1] - time[idx])

	return v_max 


def find_time_profile(time): # OUTPUTS: a time profile based on the time vector [0, T] with the accuracy of 0.01 milisecond
	time = np.array(time)
	time_new = np.arange(0, 1000*time[-1]+1, 1)/1000
	return time_new


def find_idx(time, time_new): # OUTPUTS: the index of "time-vector-elements" inside time profile
	time_idx = np.zeros(time.shape)

	counter = 0
	for i, e in enumerate(time_new):
		if e>=time[counter]:
			time_idx[counter] = i
			counter += 1

	return time_idx


def bridge_idx(velocity_profile, v_max):
	v_max_i 	= v_max[0]
	v_max_ip1 	= v_max[1]
	a = 0 
	b = velocity_profile.shape[0] - 1

	print(v_max)
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

def find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile): # OUTPUTS: the velocity profile and position profile
	time_index = find_idx(time_vector, (time_profile))
	velocity_profile = np.zeros(time_profile.shape)
	position_profile = np.zeros(time_profile.shape)


	k = 0.5 # k = T_a/T
	
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
	current_pos = position_vector[0]
	for idx, v in enumerate(velocity_profile[1:]):
		current_pos += v*time_step
		position_profile[idx+1] = current_pos

	return (velocity_profile, position_profile)


def modify_velocity_profile(time_profile, velocity_profile, v_max_vector, p0):

	position_profile = np.zeros(time_profile.shape)

	v_max_idx_vector = np.zeros(v_max_vector.shape) # index of v max in v profile
	
	# print(v_max_vector)

	counter = 0
	for idx, v_idx in enumerate(velocity_profile):
		if abs(v_idx - v_max_vector[counter]) < 10**(-10):
			v_max_idx_vector[counter] = int(idx)
			counter += 1


	for idx, v_max_idx in enumerate(v_max_vector[:-1]):
		if (v_max_vector[idx]>0) and (v_max_vector[idx+1]>0):

			a, b = bridge_idx(velocity_profile[int(v_max_idx_vector[idx]):int(v_max_idx_vector[idx+1])+3], v_max_vector[idx:idx+2])
			a += int(v_max_idx_vector[idx])
			b += int(v_max_idx_vector[idx])
			
			velocity_profile[a:b] = min(v_max_vector[idx], v_max_vector[idx+1])
		elif (v_max_vector[idx]<0) and (v_max_vector[idx+1]<0):
			a, b = bridge_idx(velocity_profile[int(v_max_idx_vector[idx]):int(v_max_idx_vector[idx+1])+3], v_max_vector[idx:idx+2])
			a += int(v_max_idx_vector[idx])
			b += int(v_max_idx_vector[idx])
			
			velocity_profile[a:b] = max(v_max_vector[idx], v_max_vector[idx+1])
		else:
			pass


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
	time_profile = find_time_profile(time_vector)																				# returns (new time profile (miliseconds), index of the time vector elements in the time profile)
	velocity_profile, position_profile = find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile)		# vector profile (trapezoidal)

	velocity_profile_new, position_profile_new = modify_velocity_profile(time_profile, np.copy(velocity_profile), v_max_vector, position_vector[0])

	plt.plot(time_profile, velocity_profile)
	plt.plot(time_profile, velocity_profile_new)
	plt.show()
	plt.plot(time_profile, position_profile)
	plt.plot(time_profile, position_profile_new)
	plt.show()

# =================================================================================================
# -------------------------------------------------------------------------------------------------
# =================================================================================================

main([0, 0.05])
main([0, 0.02, 0.05, 0.0])
main([0, 0.02, 0.03, 0.05, -0.03, 0.0])
main([0, 0.02, 0.05, 0.04, -0.02, -0.08, 0.0])
main([0, 0.02, 0.05, 0.04, -0.02, -0.06, -0.08, -0.05, 0.0])