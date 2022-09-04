
# =================================================================================================
# -- IMPORT ---------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
import matplotlib.pyplot as plt 
import time 

# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - POINT TO POINT -------------------------------------------------
# =================================================================================================

def ptp_duration(time_profile, k):
	T = time_profile[-1] - time_profile[0]
	T_a = T*k 
	return T, T_a


def acc(v_max, T_a):
	return v_max/T_a


def ptp_p_profile(time_profile, v_max, a, p_i, T, T_a):
	position_profile = np.zeros(time_profile.shape)

	n_a = 0 
	for idx, t_idx in enumerate(time_profile):
		if t_idx <= T_a:
			n_a = idx 
		else:
			break

	position_profile[0: n_a+1] 		= a/2*time_profile[0: n_a+1]**2 + p_i
	position_profile[n_a+1: -n_a-1]  	= a*T_a**2/2 + v_max*(time_profile[n_a+1: -n_a-1] - T_a) + p_i
	position_profile[-n_a-1:] 			= a*T_a**2/2 + v_max*(T - 2*T_a) + a/2*(T - T_a)**2 - a*T*(T - T_a) + a*T*time_profile[-n_a-1:] - a*time_profile[-n_a-1:]**2/2 + p_i

	return position_profile


def ptp_v_profile(time_profile, v_max, a, T, T_a):
	velocity_profile = np.zeros(time_profile.shape)

	n_a = 0 
	for idx, t_idx in enumerate(time_profile):
		if t_idx <= T_a:
			n_a = idx 
		else:
			break

	velocity_profile[0: n_a+1] 		= a*time_profile[0: n_a+1] 
	velocity_profile[n_a+1: -n_a-1]  	= v_max
	velocity_profile[-n_a-1:] 			= a*(T - time_profile[-n_a-1:])
	return velocity_profile

# =================================================================================================
# -- POINT TO POINT -------------------------------------------------------------------------------
# =================================================================================================

def ptp_velocity_profile(time_profile, position, v_max):

	p_i = position[0] 	# initial position 
	p_f = position[1] 	# final position 
	k = 0.5 			# T_a/T

	time_profile = time_profile[:] - time_profile[0]

	T, T_a = ptp_duration(time_profile, k)		# duration of point to point movement from start to finish 
	a = acc(v_max, T_a)						# returns the measure of acceleration
	# print(T, T_a, a)

	position_profile = ptp_p_profile(time_profile, v_max, a, p_i, T, T_a)		# returns the position profile
	velocity_profile = ptp_v_profile(time_profile, v_max, a, T, T_a)			# returns velocity profile

	return (velocity_profile, position_profile)

# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - MULTI POINT ----------------------------------------------------
# =================================================================================================

def find_time_vector(position):
	position = np.array(position)
	time_vector = np.linspace(0, 1, position.shape[0])
	return time_vector


def find_v_max(time, position):
	time = np.array(time)
	position = np.array(position)
	v_max = np.zeros(time.shape)

	for idx, time_idx in enumerate(time[:-1]):
		v_max[idx] = 2*(position[idx+1] - position[idx])/(time[idx+1] - time[idx])

	return v_max 


def find_time_profile(time):
	time = np.array(time)
	time_new = np.arange(0, 1000*time[-1]+1, 1)/1000
	return time_new


def find_idx(time, time_new):
	time_old_idx = np.zeros(time.shape)

	counter = 0 
	for idx, time_new_idx in enumerate(time_new):
		if abs(time_new_idx - time[counter]) < 10**(-5):
			time_old_idx[counter] = int(idx)
			counter += 1

	return time_old_idx


# =================================================================================================
# -- MULTI POINT ----------------------------------------------------------------------------------
# =================================================================================================

def find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile):

	time_index = find_idx(time_vector, time_profile)
	velocity_profile = np.zeros(time_profile.shape)
	position_profile = np.zeros(time_profile.shape)

	counter = 0
	for idx, time_idx in enumerate(time_vector[:-1]):
		temp_position 	= position_vector[idx:idx+2]										# initial and final position in time interval number i (encoded as idx)
		temp_time 		= time_vector[idx: idx+2]											# initial and final time in time interval number i (encoded as idx)
		temp_time_new 	= time_profile[int(time_index[idx]): int(time_index[idx+1])+1] 		# time profile in time interval number i (encoded as idx)

		temp_velocity, temp_position = ptp_velocity_profile(temp_time_new, temp_position, v_max_vector[idx])	# calculating the velocity profile in time interval number i (encoded as idx)

		velocity_profile[counter:counter+temp_velocity.shape[0]] = temp_velocity
		position_profile[counter:counter+temp_velocity.shape[0]] = temp_position

		counter += temp_velocity.shape[0]-1

	return (velocity_profile, position_profile)


def modify_velocity_profile(time_profile, velocity_profile_new, v_max_vector):
	
	v_max_idx_vector = np.zeros(v_max_vector.shape) # index of v max in v profile

	counter = 0
	for idx, v_idx in enumerate(velocity_profile_new):
		if abs(v_idx - v_max_vector[counter]) < 10**(-5):
			v_max_idx_vector[counter] = int(idx)
			counter += 1

	for idx, v_max_idx in enumerate(v_max_vector[:-1]):
		if (v_max_vector[idx]>0) and (v_max_vector[idx+1]>0):
			velocity_profile_new[int(v_max_idx_vector[idx]):int(v_max_idx_vector[idx+1]+1)] = min(v_max_vector[idx], v_max_vector[idx+1])
		elif (v_max_vector[idx]<0) and (v_max_vector[idx+1]<0):
			velocity_profile_new[int(v_max_idx_vector[idx]):int(v_max_idx_vector[idx+1]+1)] = max(v_max_vector[idx], v_max_vector[idx+1])
		else:
			pass

	return velocity_profile_new

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

def main():

	# constants
	position_vector = [0, 0.02, 0.05, 0.04, 0.01, 0.0]		# the position we want EE to cross

	# the vectors 
	time_vector = find_time_vector(position_vector)			# define the time vector (equally distanced times)
	v_max_vector = find_v_max(time_vector, position_vector)	# calculate v_max in each time interval

	# porfiles
	time_profile = find_time_profile(time_vector) 															# returns (new time profile (miliseconds), index of the time vector elements in the time profile)
	velocity_profile, position_profile = find_velocity_profile(time_vector, position_vector, v_max_vector, time_profile)		# vector profile (trapezoidal)

	velocity_profile_new = modify_velocity_profile(time_profile, velocity_profile, v_max_vector)


# =================================================================================================
# -------------------------------------------------------------------------------------------------
# =================================================================================================

main()