
# =================================================================================================
# -- IMPORT ---------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
import matplotlib.pyplot as plt 
import time 

# =================================================================================================
# -- MOTION KINEMATICS ----------------------------------------------------------------------------
# =================================================================================================

def ptp_duration(p_i, p_f, v_max, k):
	T = (p_f - p_i)/(v_max*(1 - k))
	return T

def acc_duration(T, k):
	T_a = T*k 
	return T_a

def acc(v_max, T_a):
	a = v_max/T_a
	return a 

def t_profile(T, n): 
	t = np.linspace(0, T, n)
	return t 

def p_profile(t, v_max, a, p_i, T, T_a):
	p_profile = np.zeros(t.shape)

	n_a = 0 
	for idx, t_idx in enumerate(t):
		if t_idx <= T_a:
			n_a = idx 
		else:
			break

	p_profile[0: n_a+1] 		= a/2*t[0: n_a+1]**2 + p_i
	p_profile[n_a+1: -n_a-1]  	= a*T_a**2/2 + v_max*(t[n_a+1: -n_a-1] - T_a) + p_i
	p_profile[-n_a-1:] 			= a*T_a**2/2 + v_max*(T - 2*T_a) + a/2*(T - T_a)**2 - a*T*(T - T_a) + a*T*t[-n_a-1:] - a*t[-n_a-1:]**2/2 + p_i

	return p_profile

def v_profile(t, v_max, a, T, T_a):
	v_profile = np.zeros(t.shape)

	n_a = 0 
	for idx, t_idx in enumerate(t):
		if t_idx <= T_a:
			n_a = idx 
		else:
			break

	v_profile[0: n_a+1] 		= a*t[0: n_a+1] 
	v_profile[n_a+1: -n_a-1]  	= v_max
	v_profile[-n_a-1:] 			= a*(T - t[-n_a-1:])
	return v_profile

# =================================================================================================
# -- MAIN ---------------------------------------------------------------------------------------
# =================================================================================================

def main():
	p_i = 1 	# initial position 
	p_f = 5 	# final position 
	v_max = 1 	# maximum velocity
	k = 0.3 	# T_a/T
	n = 50

	T = ptp_duration(p_i, p_f, v_max, k)	# duration of point to point movement from start to finish 
	T_a = acc_duration(T, k)				# duration of acceleration and deceleration 
	a = acc(v_max, T_a)						# returns the measure of acceleration 				

	time_profile = t_profile(T, n)							# returns time profile 
	position_profile = p_profile(time_profile, v_max, a, p_i, T, T_a)
	velocity_profile = v_profile(time_profile, v_max, a, T, T_a)			# returns velocity profile

	plt.plot(time_profile, velocity_profile)
	plt.show()

	plt.plot(time_profile, position_profile)
	plt.show()

# =================================================================================================
# ---------------------------------------------------------------------------------------------
# =================================================================================================

main()