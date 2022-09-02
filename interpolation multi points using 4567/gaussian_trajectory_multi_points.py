
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

def find_s_coefficients(tau, velocity, position):
    velocity = np.array(velocity)
    position = np.array(position)
    tau = np.array(tau)

    v_i = velocity[0]
    v_i_plus_1 = velocity[1]

    p_i = position[0]
    p_i_plus_1 = position[1]

    a = v_i_plus_1 + v_i + 2*p_i - 2*p_i_plus_1
    b = 3*p_i_plus_1 - 3*p_i - 2*v_i - v_i_plus_1
    c = v_i
    d = p_i

    s = np.array([a, b, c, d])
    s.reshape((1, 4))

    return s

# =================================================================================================

def v_gaussian(position_vector, v_max, center=0.5, deviation=0.18):

    position_vector = np.array(position_vector)
    t = np.linspace(0, 1, position_vector.shape[0]) # assigning t vector with equal time intervals

    a = v_max   # maximum velocity is the same as maximum height
    b = center  # because time is normalized the center must be at 0.5
    c = deviation   # standard deviation

    v_profile = a*np.exp(-(t-b)**2/(2*c**2))

    return (t, v_profile)

# =================================================================================================

def v_profile(t, v_vector, position_vector, n_segment=10):
    t =  np.array(t) 
    v_vector = np.array(v_vector)
    n = t.shape[0] - 1

    t_new = np.linspace(0, 1, n_segment*n + 1)
    
    # intializing time, position and velocity profile (and s coefficient matrix(a, b, c, d)_i)
    v_profile = np.zeros(t_new.shape)
    position_profile = np.zeros(t_new.shape)
    a_profile = np.zeros(t_new.shape)
    s = np.zeros((t.shape[0], 4)) 


    for idx, t_idx in enumerate(t[:-1]):

        tau = (t_new[idx*n_segment: idx*n_segment + n_segment] - t[idx])/(t[idx+1] - t[idx])
        s[idx, :] = find_s_coefficients(tau, v_vector[idx:idx+2], position_vector[idx:idx+2])
        
        # assigning coefficients of s_i
        a = s[idx, 0]
        b = s[idx, 1]
        c = s[idx, 2]
        d = s[idx, 3]
        position_profile[idx*n_segment: idx*n_segment + n_segment] = a*tau**3 + b*tau**2 + c*tau**1 + d

        v_profile[idx*n_segment: idx*n_segment + n_segment] = 3*a*tau**2 + 2*b*tau**1 + c
        a_profile[idx*n_segment: idx*n_segment + n_segment] = 6*a*tau**1 + 2*b
    

    return (t_new, position_profile, v_profile, a_profile)

# =================================================================================================
# -- MAIN --------------------------------------------------------------------------------------
# =================================================================================================

v_max = 0.05
position_vector = [0, 0.01, 0.03, 0.07, 0.09, 0.10, 0.15, 0.18]
t, v_vector = v_gaussian(position_vector, v_max)

t_new, position_profile, v_profile, a_profile = v_profile(t, v_vector, position_vector)


print("--- %s seconds ---" % (time.time() - start_time))

# plt.plot(t, position_vector)
# plt.plot(t_new, position_profile)
# plt.show()


plt.plot(t, v_vector)
plt.plot(t_new, v_profile)
plt.show()

# plt.plot(t_new, a_profile)
# plt.show()





