import gl.gl as gl
import at.at as atm
import vsp.vsp as vsp
import vs.vs as vs
import alg.alg as alg
import sys, re, math, os, copy, commands

class myclass:
	def __init__(self, control_file):
		for i in open(control_file, "r").readlines():
			s = i.split()
			if len(s) >= 2:
				if s[0] == "poscar":
					self.poscar = s[1]
				elif s[0] == "out_poscar":
					self.out_poscar = s[1]
				elif s[0] == "vc":
					self.vc = float(s[1])
				elif s[0] == "ssize":
					self.ssize = [int(s[1]), int(s[2]), int(s[3])]
				elif s[0] == "tangent":
					self.theta = -math.atan2(float(s[3])*float(s[4]), float(s[1])*float(s[2]))
				
	def input_poscar(self):
		self.at = vsp.input_poscar(self.poscar)
	
	def frac_to_cart(self):
		self.at = atm.frac_to_cart(self.at, "vasp")
	
	def supercell(self):
		self.at = atm.supercell(self.at, self.ssize, "vasp")
	
	def cart_to_frac(self):
		self.at = atm.cart_to_frac(self.at, "vasp")
	
	def add_vc(self):
		self.at = atm.frac_to_cart(self.at, "vasp")
		
		self.at.lt[0] += self.vc
		self.at.lv = atm.ltlv(self.at.lt, "vasp")
		
		self.at = atm.cart_to_frac(self.at, "vasp")
		
	def output_poscar(self):
		for i in range(self.at.atn):
			if self.at.e[i] == "Y":
				self.at.e[i] = "Al"
		
		self.at.el = ["Al"]
		vsp.output_poscar(self.at, self.out_poscar)
		
	def output_xyz(self, fn):
		self.at = atm.frac_to_cart(self.at, "vasp")
		vs.output_xyz(self.at, fn)
	
	def make_supercell(self):
		e, x, y, z = [], [], [], []
		xl, yl = [-10, 10], [-10, 10]
		
		for i in range(xl[0], xl[1]):
			for j in range(yl[0], yl[1]):
				for l in range(self.at.atn):
					e.append(self.at.e[l])
					z.append(self.at.z[l])
					x.append(self.at.x[l] + i)
					y.append(self.at.y[l] + j)
		
		self.at.e, self.at.x = list(e), list(x)
		self.at.y, self.at.z = list(y), list(z)
		self.at.atn = self.at.atn * (xl[1]-xl[0]) * (yl[1]-yl[0])
	
	def rotation_xy(self):
		delta = 1e-4
		s = math.sin(self.theta)
		c = math.cos(self.theta)
		
		for i in range(self.at.atn):
			xx = self.at.x[i]*c-self.at.y[i]*s
			yy = self.at.x[i]*s+self.at.y[i]*c
			self.at.x[i], self.at.y[i] = xx, yy
			
			if -delta <= self.at.x[i] <= delta: self.at.x[i] = 0.0
			if -delta <= self.at.y[i] <= delta: self.at.y[i] = 0.0
		
	def find_boundary_atom(self):
		delta = 1e-4
		temp_lt = set([])
		ref_element = "Y"
		
		for i in range(self.at.atn):
			if self.at.e[i] == ref_element:
				if self.at.y[i] >= delta:
					if -delta <= self.at.x[i] <= delta:
						if -delta <= self.at.z[i] <= delta:
							temp_lt.add(self.at.y[i])
						
		temp_lt = list(temp_lt)
		
		if len(temp_lt) == 1:
			self.at.lt[1] = temp_lt[0]
		elif len(temp_lt) == 0:
			sys.exit("No boundary atom is found in the y direction.")
		elif len(temp_lt) >= 2:
			print "Several boundary atoms are found in the x direction."
			print "The smallest value is set to self.at.lt[1]."
			self.at.lt[1] = min(temp_lt)
		print "self.at.lt[1] = %f" % (self.at.lt[1])
		
		temp_lt = set([])
		for i in range(self.at.atn):
			if self.at.e[i] == ref_element:
				if self.at.x[i] >= delta:
					if -delta <= self.at.y[i] <= delta:
						if -delta <= self.at.z[i] <= delta:
							temp_lt.add(self.at.x[i])
						
		temp_lt = list(temp_lt)		
		
		if len(temp_lt) == 1:
			self.at.lt[0] = temp_lt[0]
		elif len(temp_lt) == 0:
			sys.exit("No boundary atom is found in the x direction.")
		elif len(temp_lt) >= 2:
			print "Several boundary atoms are found in the x direction."
			print "The smallest value is set to self.at.lt[0]."
			self.at.lt[0] = min(temp_lt)
		print "self.at.lt[0] = %f" % (self.at.lt[0])

	def cut_cell(self):
		delta = 1e-4
		e, x, y, z = [], [], [], []
		
		self.at.lv = atm.ltlv(self.at.lt, "vasp")
		self.at = atm.cart_to_frac(self.at, "vasp")

		for i in range(self.at.atn):
			if -delta <= self.at.x[i] <= 1.0-delta:
				if -delta <= self.at.y[i] <= 1.0-delta:
					if -delta <= self.at.z[i] <= 1.0-delta:
						e.append(self.at.e[i])
						x.append(self.at.x[i])
						y.append(self.at.y[i])
						z.append(self.at.z[i])
		
		self.at.atn = len(e)
		self.at.e, self.at.x = list(e), list(x)
		self.at.y, self.at.z = list(y), list(z)
	
	def rotation(self):
		for i in range(self.at.atn):
			self.at.y[i] = 1. - self.at.y[i]
			self.at.z[i] = 1. - self.at.z[i]
	
	def wrap_atom(self):
		for i in range(self.at.atn):
			if self.at.x[i] >= 0.9999: self.at.x[i] -= 1.
			elif self.at.x[i] < 0.: self.at.x[i] += 1.
			
			if self.at.y[i] >= 0.9999: self.at.y[i] -= 1.
			elif self.at.y[i] < 0.: self.at.y[i] += 1.
			
			if self.at.z[i] >= 0.9999: self.at.z[i] -= 1.
			elif self.at.z[i] < 0.: self.at.z[i] += 1.
	
	def sw_yz(self):
		self.at.y, self.at.z = list(self.at.z), list(self.at.y)
		self.at.lt[1], self.at.lt[2] = self.at.lt[2], self.at.lt[1]
		self.at.lv = atm.ltlv(self.at.lt, "vasp")
	
	def sw_xz(self):
		self.at.x, self.at.z = list(self.at.z), list(self.at.x)
		self.at.lt[0], self.at.lt[2] = self.at.lt[2], self.at.lt[0]
		self.at.lv = atm.ltlv(self.at.lt, "vasp")
	
################################################################################
control_file = "control_110.in"

ob = myclass(control_file)

ob.input_poscar()
#ob.adjust_position()
ob.make_supercell()
ob.frac_to_cart()
ob.rotation_xy()
ob.find_boundary_atom()
ob.cut_cell()

ob.supercell()
ob.add_vc()
ob.sw_yz()
ob.output_poscar()
