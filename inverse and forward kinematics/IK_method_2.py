class Jacobian:
	def __init__(self, EE_position, active_rod=0.2, passive_rod=0.46, base_radius=0.3464101615, EE_radius=0.2563435195, alpha=[0, 120, 240]):
		# i made a mistake during naming of the variables, base and EE "radius" are not exactly radiuses, they are sides of a triangle ...
		# ... made by EE joints (the joints are the middle of the side of the traingle)
		
		# initializing the basic geometry and the given data
		self.alpha = np.array(alpha)							# alpha angles
		self.EE_position_global = np.array(EE_position)			# end effoctor position (x_e, y_e, z_e) with respect to alpha1 = 0								
		self.active_rod = active_rod							# length of the active rod (the upper rod or r_f)
		self.passive_rod = passive_rod							# length of the passive rod (the lower rod or r_e)
		self.EE_radius = EE_radius								# the radius of the end effoctor e 
		self.base_radius = base_radius							# the radius of the base or f

	def get_theta_ij(self): 
		# this also calculate inverse kinematic

		# assigning constants
		alpha = math.pi/180*self.alpha
		R = self.base_radius*(3**0.5/6)
		r = self.EE_radius*(3**0.5/6)
		a = self.active_rod
		b = self.passive_rod

		# assigning position of the EE 
		px = self.EE_position_global[0]
		py = self.EE_position_global[1]
		pz = self.EE_position_global[2]

		# initializing theta 1, 2, 3
		theta_1 = np.zeros((3))
		theta_2 = np.zeros((3))
		theta_3 = np.zeros((3))

		# calculating theta 1, 2, 3
		for i in [0, 1, 2]:
			theta_3[i] = math.acos((px*sin(alpha[i]) + py*cos(alpha[i]))/b)

			A = px*cos(alpha[i]) - py*sin(alpha[i]) - R + r
			B = pz
			M = (A**2 + B**2 + a**2 - (b*sin(theta_3[i]))**2)/(2*a)
			t = (B + (B**2 - M**2 + A**2)**0.5)/(M + A)
			theta_1[i] = 2*math.atan(t)
			theta_2[i] = math.asin((pz - a*sin(theta_1[i]))/(b*sin(theta_3[i]))) - theta_1[i]
		
		# now we have all of the theta 1, 2 and 3 
		self.theta_1 = theta_1
		self.theta_2 = theta_2
		self.theta_3 = theta_3
		return theta_1
	
	def get_jacobian_matrix(self): 
		
		# initializing J_ij 
		jx = np.zeros((3))
		jy = np.zeros((3))
		jz = np.zeros((3))
		J_P = np.zeros((3, 3))
		J_theta = np.zeros((3, 3))

		for i in [0, 1, 2]:
			jx[i] =  sin(self.theta_3[i])*cos(self.theta_2[i] + self.theta_1[i])*cos(self.alpha[i]) + cos(self.theta_3[i])*sin(self.alpha[i])
			jy[i] = -sin(self.theta_3[i])*cos(self.theta_2[i] + self.theta_1[i])*sin(self.alpha[i]) + cos(self.theta_3[i])*cos(self.alpha[i])
			jz[i] =  sin(self.theta_3[i])*sin(self.theta_2[i] + self.theta_1[i])
			J_P[i, :] = [jx[i], jy[i], jz[i]]
			J_theta[i, i] = sin(self.theta_2[i])*sin(self.theta_3[i])
		
		return (J_P, J_theta)
