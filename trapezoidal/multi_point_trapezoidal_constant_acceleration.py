
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



# =================================================================================================
# -- POINT TO POINT -------------------------------------------------------------------------------
# =================================================================================================



# =================================================================================================
# -- MOTION KINEMATICS FUNCTIONS - MULTI POINT ----------------------------------------------------
# =================================================================================================



# =================================================================================================
# -- MULTI POINT ----------------------------------------------------------------------------------
# =================================================================================================


def find_v_max(time, position):	# OUTPUTS: the v_max in time intervals [t_i, t_(i+1)]
	time = np.array(time)
	position = np.array(position)
	v_max = np.zeros(time.shape[0])

	for idx, time_idx in enumerate(time[:-1]):
		v_max[idx] = 2*(position[idx+1] - position[idx])/(time[idx+1] - time[idx])

	return v_max 


# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

def main(position_vector):

	a = 0.2

	# the vectors 
	v_max_vector = find_v_max(time_vector, position_vector)	# calculate v_max in each time interval
	time_vector = find_time_vector(position_vector)			# define the time vector (equally distanced times)
	

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