
# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import matplotlib.pyplot as plt 
import time 
import math 

# =================================================================================================
# -- FORWARD KINEMATICS ---------------------------------------------------------------------------
# =================================================================================================
# INPUT: actuator angles (unit degrees) 
# OUTPUT: position of EE (unti meters)

def forward_kinematics(Theta, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):

		alpha = np.array(alpha)								# alpha angles
		Theta = np.array(Theta)								# actuator angles in degrees

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
# -- TEST -----------------------------------------------------------------------------------------
# =================================================================================================

print(forward_kinematics([-6.11174185, -6.11174185, -6.11174185]))
print(forward_kinematics([ 41.78484676, -18.12320943,  10.51376895]))
