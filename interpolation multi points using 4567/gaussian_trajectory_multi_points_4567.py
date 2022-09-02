
# =================================================================================================
# -- IMPORT ---------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
import matplotlib.pyplot as plt 
import time 

start_time = time.time()

# =================================================================================================
# -- MOTION PLANNING ------------------------------------------------------------------------------
# =================================================================================================

def time_vector(position_vector): # done
	position_vector = np.array(position_vector)
	t = np.linspace(0, 1, position_vector.shape[0])   # assigning t vector with equal time intervals

	return t
	

def v_gaussian(t, v_max, center, deviation): # done

	t = np.array(t)
	
	a = v_max   # maximum velocity is the same as maximum height
	b = center  # because time is normalized the center must be at 0.5
	c = deviation   # standard deviation

	v_vector = a*np.exp(-(t-b)**2/(2*c**2))

	return v_vector

def a_gaussian(t, v_max, center, deviation): # done

	t = np.array(t)
	
	a = v_max   # maximum velocity is the same as maximum height
	b = center  # because time is normalized the center must be at 0.5
	c = deviation   # standard deviation

	M = np.exp(-0.5*((t-b)/c)**2)
	a_vector = (b-t)/c**2*a*M

	return a_vector

def j_gaussian(t, v_max, center, deviation): # done

	t = np.array(t)
	
	a = v_max   # maximum velocity is the same as maximum height
	b = center  # because time is normalized the center must be at 0.5
	c = deviation   # standard deviation

	M = np.exp(-0.5*((t-b)/c)**2)
	j_vector = a/c**2*(((b - t)**2/c**2) - 1)*M

	return j_vector

# =================================================================================================
# -- INTERPOLATION --------------------------------------------------------------------------------
# =================================================================================================

def inter_coeff_matrix(t, p, v, a, j):
	# converting to numpy 
	t = np.array(t)
	p = np.array(p)
	v = np.array(v)
	a = np.array(a)
	j = np.array(j)


	# initializing s 
	s = np.zeros((t.shape[0], 8))

	for idx, t_idx in enumerate(t[:-1]):

		# p_i_1 means p_(i+1)
		p_i = p[idx]
		p_i_1 = p[idx+1]
		v_i = v[idx]
		v_i_1 = v[idx+1]
		print(idx)
		print(a)
		a_i = a[idx]
		a_i_1 = a[idx+1]
		j_i = j[idx]
		j_i_1 = j[idx+1]

		a_coeff = 4*a_i - 2*a_i_1 + j_i + j_i_1/6 + 20*p_i - 20*p_i_1 + 10*v_i + 10*v_i_1
		b = 13/2*a_i_1 - 15*a_i - 4*j_i - j_i_1/2 - 70*p_i + 70*p_i_1 - 36*v_i - 34*v_i_1
		c = 20*a_i - 7*a_i_1 + 6*j_i + j_i_1/2 + 84*p_i - 84*p_i_1 + 45*v_i + 39*v_i_1
		d = 5/2*a_i_1 - 10*a_i - 4*j_i - j_i_1/6 - 35*p_i + 35*p_i_1 - 20*v_i - 15*v_i_1
		e = j_i
		f = a_i
		g = v_i
		h = p_i

		s[idx, :] = np.array([a_coeff, b, c, d, e, f, g, h])

	return s 

# =================================================================================================
# -- PROFILES -------------------------------------------------------------------------------------
# =================================================================================================

def t_profile(t, n_segment):
	t =  np.array(t) 
	n = t.shape[0] - 1
	t_new = np.linspace(0, 1, n_segment*n + 1)

	return  t_new

def p_profile(t, t_new, v, s):
	# converting to numpy array
	t = np.array(t)
	t_new = np.array(t)
	v = np.array(v)
	s = np.array(v)

	# initializing the velocity profile
	p_profile = np.zeros(t_new.shape)

	return p_profile

def v_profile(t, t_new, v, s, n_segment):
	# converting to numpy array
	t = np.array(t)
	t_new = np.array(t_new)
	v = np.array(v)
	s = np.array(s)

	# initializing the velocity profile
	v_profile = np.zeros(t_new.shape)

	for idx, t_idx in enumerate(t[:-1]):

		tau = (t_new[idx*n_segment: idx*n_segment + n_segment] - t[idx])/(t[idx+1] - t[idx])

		# assigning coefficients of s_i
		a = s[idx, 0]
		b = s[idx, 1]
		c = s[idx, 2]
		d = s[idx, 3]
		e = s[idx, 4]
		f = s[idx, 5]
		g = s[idx, 6]
		h = s[idx, 7]

		v_profile[idx*n_segment: idx*n_segment + n_segment] = 7*a*tau**6 + 6*b*tau**5 + 5*c*tau**4 + 4*d*tau**3 + 3*e*tau**2 + 2*f*tau**1 + 1*g

	return v_profile 

def a_profile(t, t_new, v, s):


	return a_profile 

def j_profile(t, t_new, v, s):


	return j_profile 

# =================================================================================================
# -- MAIN -----------------------------------------------------------------------------------------
# =================================================================================================

def main():

	v_max = 1
	n_segment = 1000
	position_vector = [0, 0.01, 0.03, 0.07, 0.09, 0.10, 0.15, 0.18]
	t_vector = time_vector(position_vector) 												# normalized time vector with equal time intervals (discrete values)

	center_of_curve = 0.5
	standard_deviation = 0.1

	# =========================== discreting vectors of the data ===========================
	velocity_vector = v_gaussian(t_vector, v_max, center_of_curve, standard_deviation) 		# returns the overall velocity vector (discrete values)
	acceleration_vector = a_gaussian(t_vector, v_max, center_of_curve, standard_deviation) 	# returns the overall acceleration vector (discrete values)
	jerk_vector = j_gaussian(t_vector, v_max, center_of_curve, standard_deviation)			# returns the overall jerk vector (discrete values)


	# =========================== interpolating the vectors and gaining the profiles ===========================
	t_new = t_profile(t_vector, n_segment)																	# returns time profile

	s = inter_coeff_matrix(t_vector, position_vector, velocity_vector, acceleration_vector, jerk_vector)

	# position_profile = p_profile(t_vector, t_new, velocity_vector, s, n_segment)
	velocity_profile = v_profile(t_vector, t_new, acceleration_vector, s, n_segment)
	# acceleration_profile = a_profile(t_vector, t_new, jerk_vector, s, n_segment)
	# jerk_profile = j_profile(t_vector, t_new, jerk_vector, s, n_segment)


	# =========================== ploting the gethered data ===========================
	# plt.plot(t_vector, position_vector)				# ploting the discrete position 
	# plt.plot(t_new, position_profile)				# ploting the position profile 
	# plt.show()										# comparing the position vector with its interpolation


	plt.plot(t_vector, velocity_vector)				# ploting the discrete velocity 
	plt.plot(t_new, velocity_profile)				# ploting the velocity profile 
	plt.show()										# comparing the velocity vector with its interpolation

	# plt.plot(t_vector, acceleration_vector)			# ploting the discrete acceleration 
	# plt.plot(t_new, acceleration_profile)			# ploting the acceleration profile 
	# plt.show()										# comparing the acceleration vector with its interpolation

	# plt.plot(t_vector, jerk_vector)					# ploting the discrete jerk 
	# plt.plot(t_new, jerk_profile)					# ploting the jerk profile 
	# plt.show()										# comparing the jerk vector with its interpolation

# =================================================================================================
# -------------------------------------------------------------------------------------------------
# =================================================================================================

main()