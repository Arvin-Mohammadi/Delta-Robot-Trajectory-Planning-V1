# attempt #1 
# trajectory planning, point to point movement with 5rd order polynomial 
# 2 fundamental questions: 
# 1. given a specific location in the real world, what values should my robot's joint be ...
# ... set to in order to get the EE there? (inverse kinematics)
# 2. given the setting of my joints, where is my EE in real world coordinates? (forward kinematics)

# =================================================================================================
# -- imports --------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
from sympy import *

# =================================================================================================
# -- inverse kinematic class ----------------------------------------------------------------------
# =================================================================================================

# here we have the EE position and the lengths of rods and basic geometry
# our goal is to find theta1, theta2 and theta3 

class InverseKinematics:
	def __init__(self, EE_position, active_rod, passive_rod, base_radius, EE_raidus):
		# initializing the basic geometry and the given data

		self.EE_position = np.array(EE_position)	# end effoctor position (x_e, y_e, z_e)
		self.active_rod = active_rod		# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod		# length of the passive rod (the lower rod or r_e)
		self.EE_raidus = EE_raidus			# the radius of the end effoctor e
		self.base_radius = base_radius		# the radius of the base or f

	def get_positions(self):
		# here positions of the the 3 important points are found E1, E1_prime, F1
		
		self.E1_position = self.EE_position + np.array([0, -self.EE_raidus/3, 0]) 
		self.E1_prime_position = np.array([0, self.E1_position[1], self.E1_position[2]])
		self.F1_position = np.array([0, - self.base_radius/3, 0])

		# and from those 3 point we find J1 
		# (intersection of 2 circles, with centers of F1 and E1_prime and ridus of r_f and (r_e**2 - E1_x**2)**0.5)
		# equation is :
		# (r_f**2 - (y - y_F)**2)**0.5 + z_F == (r_e**2 - x_E1**2 - (y - y_E1prime)**2)**0.5 + z_E1prime
		r_f = float(self.active_rod)
		r_e = float(self.passive_rod)
		y_F = float(self.F1_position[1])
		z_F = float(self.F1_position[2])
		x_E1 = float(self.E1_position[0])
		y_E1prime = float(self.E1_prime_position[1])
		z_E1prime = float(self.E1_prime_position[2])

		y = symbols('y')
		z = symbols('z')

		equation1 = Eq((y - y_F)**2 + (z - z_F)**2 - r_f**2, 0)
		equation2 = Eq((y - y_E1prime)**2 + (z - z_E1prime)**2 - (r_e**2 - x_E1**2), 0)
		solution = solve((equation1, equation2), (y, z))

		[y, z] = solution[0]
		self.J1_position = [0, y, z]

	def get_theta(self):
		z_J1 = self.J1_position[2]
		y_J1 = self.J1_position[1]
		y_F1 = self.F1_position[1]

		self.theta1 = math.atan(z_J1/(y_F1 - y_J1))
		print(self.theta1)
		return self.theta1

# =================================================================================================
# -- main -----------------------------------------------------------------------------------------
# =================================================================================================


Kinematics = InverseKinematics([0, 0, -0.5], 0.1, 0.4, 0.3, 0.3)
Kinematics.get_positions()
Kinematics.get_theta()


