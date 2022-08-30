# B-Spline

# =================================================================================================
# -- IMPORTS --------------------------------------------------------------------------------------
# =================================================================================================

import numpy as np 
import math 
import matplotlib.pyplot as plt 
from numpy import matmul 
from numpy.linalg import inv
np.set_printoptions(precision = 2)

# =================================================================================================
# -- B SPLINE --------------------------------------------------------------------------------------
# =================================================================================================

# this function only works for p=4 

class BSpline: 
	def __init__(self, q, t, T=0.1, vi=0, vf=0, ai=0, af=0): 
		self.q = np.array(q)	# the points coordinates 
		self.t = np.array(t)	# the time instants array 
		self.n = self.q.shape[0] - 1

	def get_u_vector(self, t, p=4):
		t = np.array(t)
		n = t.shape[0] - 1

		if p%2 == 0: # p even
			n_knot = n + 2*p + 1
			u = np.zeros((n_knot+1))
			u[0:p+1] = t[0]
			u[p+1:n+p+1] = (t[0:n] + t[1:n+1])/2
			u[n+p+1:n+2*p+2] = t[n]
		elif p%2 != 0: # p odd 
			n_knot = n + 2*p
			u = np.zeros((n_knot+1))
			u[0:p+1] = t[0]
			u[p+1:n+p] = t[1:n]
			u[n+p:n+2*p+1] = t[n]
		else:
			print("p needs to be a positive integer")

		return u

	def which_span(self, u, u_instant, p=4):
		u = np.array(u)
		n_knot = u.shape[0] - 1
		high = n_knot - p 
		low = p + 1 

		if u_instant == u[high]:
			mid = high -1
		else: 
			for j, u_j in enumerate(u):
				if (u_instant >= u[j]) and (u_instant < u[j+1]):
					mid = j

		return mid

	def get_basis_function(self, u, u_instant, i, p=4): 
		u = np.array(u)
		n_knot = u.shape[0] - 1
		B = np.zeros((p+1))
		DL = np.zeros((p+1))
		DR = np.zeros((p+1))

		B[0] = 1
		for j in range(1, p+1):
			DL[j] = u_instant - u[i+1-j]
			DR[j] = u[i+j] - u_instant
			acc = 0 
			for r in range(j):
				temp = B[r]/(DR[r+1] + DL[j-r])
				B[r] = acc + DR[r+1]*temp
				acc = DL[j-r]*temp
			B[j] = acc

		return B

	def get_basis_function_derivative(self, u, u_instant, i, n=3, p=4):
		Du = np.zeros((p+1, p+1))
		Ders = np.zeros((n+1, p+1))
		DL = np.zeros((p+1))
		DR = np.zeros((p+1))
		a = np.zeros((n+1, n+1))

		Du[0, 0] = 1
		for j in range(1, p+1):
			DL[j] = u_instant - u[i+1-j]
			DR[j] = u[i+j] - u_instant
			acc = 0 
			for r in range(j):
				Du[j, r] = DR[r+1] + DL[j-r]

				temp = Du[r, j-1]/Du[j, r]

				Du[r, j] = acc + DR[r+1]*temp
				acc = DL[j-r]*temp 
			Du[j, j] = acc 

		for j in range(p+1):
			Ders[0, j] = Du[j, p]

		for r in range(p+1):
			s1 = 0 
			s2 = 1 
			a[0, 0] = 1

			for k in range(1, n+1):
				d = 0 
				rk = r - k 
				pk = p - k 
				if r>=k:
					a[s2, 0] = a[s1, 0]/Du[pk+1, rk]
					d = a[s2, 0]*Du[rk, pk]
				if rk>= -1:
					j1 = 1
				else: 
					j1 = -rk
				if (r-1 <= pk):
					j2 = k-1
				else:
					j2 = p - r 
				for j in range(j1, j2+1):
					a[s2, j] = (a[s1, j] - a[s1, j-1])/Du[pk+1, rk+j]
					d += a[s2, j]*Du[rk+j, pk]

				if r<=pk:
					a[s2, k] = -a[s1, k-1]/Du[pk+1, r]
					d += a[s2, k]*Du[r, pk]
				Ders[k, r] = d
				j = s1
				s1 = s2
				s2 = j 

		r = p 
		for k in range(1, n+1):
			for j in range(p+1):
				Ders[k, j] *= r 
			r*= (p-k)

		return Ders


	def get_A_matrix(self, u, p=4):
		n = self.n
		t = self.t
		n_knot = u.shape[0] - 1
		m = n_knot-p-1
		u = np.array(u)
		A = np.zeros((n+5, m+1))

		row = 0
		idx = 0
		while idx<t.shape[0]:

			i = self.which_span(u, t[idx])
			if idx == 0:
				B = traj.get_basis_function_derivative(u, t[idx], i, n=2, p=p)
				print(B)
				A[row:row+3, idx:idx+p+1] = B
				row += 2
			elif idx == t.shape[0]-1:
				print(i)
				B = traj.get_basis_function_derivative(u, t[idx], i, n=2, p=p)
				B = np.flipud(B)
				A[row:row+3, idx:idx+p+1] = B
			else:
				B = traj.get_basis_function(u, t[idx], i, p=p)
				A[row, idx:idx+p+1] = B
			
			idx += 1
			row += 1

		return A 

	def get_c_matrix(self, vi, ai, vf, af): 
		q = self.q
		n = self.n
		c = np.zeros((n+5, 1))
		c[0:3, 0] = [q[0], vi, ai]
		c[3:n+2, 0] = q[1:n]
		c[n+2:, 0] = [af, vf, q[-1]]

		return c 

	def get_P(self, A, c):
		return matmul(inv(A), c)

# =================================================================================================
# -- POSITION GENERATOR ---------------------------------------------------------------------------
# =================================================================================================

# in this part we try to generate a circle and get n sample points from that circle. 
# this results in having n-1 cartesian points (matrix.shape=(n-1, 3))

class PositionGenerator:

	def __init__(self, ratio, center):
		# t is in seconds
		n = 500
		pi = math.pi
		self.n = int(n)
		self.ratio = ratio		# r
		self.center = center	# xc, yc, zc
		self.gamma = np.linspace(0, 2*pi, num=self.n+1)

	def cart_position(self):
		# self.points is the positions of the n points in x, y and z directions (n*3 matrix)
		self.points = np.array([np.cos(self.gamma)*self.ratio + np.ones((self.n+1))*self.center[0], np.sin(self.gamma)*self.ratio + np.ones((self.n+1))*self.center[1], np.ones((self.n+1))*self.center[2]])
		self.points = np.transpose(self.points)
		return self.points

# =================================================================================================
# -- MAIN --------------------------------------------------------------------------------------
# =================================================================================================
# u = [0, 0, 0, 0, 1, 2, 4, 7, 7, 7, 7]
# u_instant = 7

# traj = BSpline([0], [0])
# i = traj.which_span(u, u_instant, p=3)
# # B = traj.get_basis_function(u, u_instant, i, p=3)
# # print(B)
# # print(i)
# B_ders = traj.get_basis_function_derivative(u, u_instant, i, n=3, p=3)
# # print(B_ders)
# # =================================================================================================
# t = [0, 5, 7, 8, 10, 15, 18]
# q = [3, -2, -5, 0, 6, 12, 8]
# vi = 2
# vf = -3
# ai = 0
# af = 0
# traj = BSpline(q, t)
# u = traj.get_u_vector(t)
# A = traj.get_A_matrix(u)
# c = traj.get_c_matrix(vi, ai, vf, af)
# print(A)
# print(c)
# p = traj.get_P(A, c)
# print(p)
# =================================================================================================
generator = PositionGenerator(0.3, [0, 0, 0])
q = generator.cart_position()[:, 0]
vi = 0 
ai = 0 
vf = 0 
af = 0 

c = np.zeros((q.shape[0] + 4))
