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
		# and from those we find J1
		self.E1_position = 

# =================================================================================================
# -- main -----------------------------------------------------------------------------------------
# =================================================================================================


Kinematics = InverseKinematics([-0.5, 0, 0], 0.3, 0.5, 0.3, 0.1)



