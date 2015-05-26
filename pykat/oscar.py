"""
------------------------------------------------------
Pyhton implementation of FFT propogration of beams.
Based on OSCAR
http://iopscience.iop.org/1742-6596/228/1/012021
http://www.mathworks.com/matlabcentral/fileexchange/20607-oscar

Work in progress, currently these functions are
untested!

Jerome, Andreas 15.04.2014
http://www.gwoptics.org/pykat/
------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import math
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import minimize
from scipy import fftpack

from numpy.polynomial.hermite import hermval



class cavity(object): # Class to represent the 2D wavefront distortion (or surface)

	def __init__(self, _inputMirror, _endMirror, _length, _inputField):
		
		self.I_input = _inputMirror
		self.I_end = _endMirror
		self.Length = _ength
		self.E_in = _inputField
		
		self.Resonance_phase = None
		self.Cavity_param = [200] # list of parameters for the running of the function
		# [Nb round trip to find the resonance phase]
		
		# Prepare the space:
		self.Field_circ = None
		self.Field_ref = None
		self.Field_trans = None
		self.Field_reso_guess = None
		self.Loss_RTL = None
		
		# To speed up, pre-compute some quantities:
		self.Propa_mat_cav = np.exp( 1j*( -self.E_in.k_prop * self.Length + \
					np.pi*(self.E_in.wavelength/self.E_in.ref_index) * \
					(self.E_in.grid.fft_X**2 + self.E_in.Grid.fft_Y**2)* self.Length ))
		
		# Pre-calculate the wavefront distortions
		self.WF_I_input = np.exp(-1j * self.E_in.k_prop * 2 * self.I_input.surface) * self.I_input.mask * self.I_input.r
		self.WF_I_end = np.exp(-1j * self.E_in.k_prop * 2 * self.I_end.surface) * self.I_end.mask * self.I_end.r
		
		
	def do_round_trip(self, Ein, RT_phase = 1): # suppose the beam at the input mirror moving toward the end mirror
		Ein.quick_propagation(self.Propa_mat_cav)
		Ein.quick_reflection(self.WF_I_end) 
		Ein.quick_propagation(self.Propa_mat_cav)
		Ein = Ein * RT_phase
		Ein.quick_reflection(self.WF_I_input) 
		return Ein
		
	def find_resonance(self): # Find the resonance length of the cavity 
		# First make the input beam pass through the input mirror
		Ein = self.E_in.copy()
		Ein.change_n(new_n = self.I_input.n2)
		Ein = Ein.ref_trans_mirror(Iin = self.I_input)[1] # Only return the transmitted field
		
		num_iter = self.Cavity_param[0]
		
		Field_total = Ein.copy()
		Field_total.normalise(power = 0)
		Field_circ = Ein.copy()
	
		Field_total.grid == Field_circ.grid
		
		Phase_adjust =1
		ii = 0
		
		while ii <  num_iter:
			Field_total = Field_total + Field_circ

			# Do a round trip for Field_circ
			Field_circ = self.do_round_trip(Field_circ,  RT_phase =  Phase_adjust)			
			Phase_adjust = Phase_adjust * np.exp(-1j*np.angle(field.amplitude.calculateOverlap(Field_circ, Field_total )))
			ii+= 1
		
		# We have the pseudo eigen mode of the cavity, will make one more round trip to find the round trip phase 
		Field_before_RT = Field_total.copy()
		Field_total = self.do_round_trip(Field_total)
		
		self.Resonance_phase = np.exp(-1j* np.angle(field.calculateOverlap(Field_total, Field_before_RT)))
		
		return self
	
	
	def calculateFields(self): # Find the circulating, transmitted and reflected fields (and also the round trip loss)
		
		# Check the number of iteration to reach the convergence
		Accuracy = 0.0001
		RT_loss = self.I_input.r* self.I_end.r
		# Have to solve RT_loss^num_iter < 0.5*accuracy
		num_iter = np.round( np.log(0.5*Accuracy)/(np.log(RT_loss)) )

		# Calculate the field entering the cavity
		Ein = self.E_in.copy()
		Ein.change_n(new_n = self.I_input.n2)
		E_dref, Ein = Ein.Ref_Trans_mirror(Iin = self.I_input)
	
		# Initialize the fields
		Field_total = Ein.copy()
		Field_total.normalise(power = 0)
		Field_circ = Ein.copy()
	
		ii = 0
		while ii <  num_iter:
			#print(ii,Field_total.Calculate_power())
			Field_total = Field_total + Field_circ

			# Do a round trip for Field_circ
			Field_circ = self.do_round_trip(Field_circ, RT_phase =  self.Resonance_phase)
			# Then adjust for the phase to be on resonance
			ii+= 1
	
		self.Field_circ = Field_total
	
		# Now calculate the transmitted and reflected field

		E_temp = Field_total.copy()
		E_temp.Quick_propagation(self.Propa_mat_cav)
		self.Field_trans = E_temp.Ref_Trans_mirror(Iin = self.I_end)[1]

		E_temp = Field_total.copy()
		E_temp.Quick_propagation(self.Propa_mat_cav)
		E_temp.Quick_reflection(self.WF_I_end) 
		E_temp.Quick_propagation(self.Propa_mat_cav)
		E_temp = E_temp * self.Resonance_phase
		E_temp2 = E_temp.Ref_Trans_mirror(Iin = self.I_input)[1]
		
		self.Field_ref = E_temp2 + E_dref
		
		# Calculate the round trip loss, take the circulating field and do one round trip
		# with the reflectivity of the mirrors set to 1
		E_temp = Field_total.copy()
		E_temp.normalise(power = 1)
		
		E_temp.Quick_propagation(self.Propa_mat_cav)
		E_temp = E_temp.Ref_Trans_mirror(Iin = self.I_end, Reflectivity = 1)[0]
		E_temp.Quick_propagation(self.Propa_mat_cav)
		E_temp = E_temp.Ref_Trans_mirror(Iin = self.I_input, Reflectivity = 1)[0]

		self.Loss_RTL = 1 - E_temp.power()
				
		return self
	
	def Display_results(self): # Display the results

		if self.Field_circ == None:
			print('Display_results(): the method calculateFields () must be run first')
		else:
			str = 'Power in the input beam: {0:10.4f} [W]'.format(self.E_in.power());print(str)
			str = 'Circulating power: {0:10.4f} [W]'.format(self.Field_circ.power());print(str)
			str = 'Transmitted power: {0:10.4f} [W]'.format(self.Field_trans.power());print(str)
			str = 'Reflected power: {0:10.4f} [W]'.format(self.Field_ref.power());print(str)
			str = 'Round trip losses: {0:10.4f} [ppm]'.format(self.Loss_RTL*1E6);print(str)
			
			# If someone know how to plot 2 E_field (with E_plot) next to each other
			# I am interested!
 			
	
	def Scan_cavity(self): # Scan the cavity over one FSR
		# First check how many round trip to do to have a reasonable idea of the circulating power
		Accuracy = 0.01
		RT_loss = self.I_input.r* self.I_end.r
		num_iter = np.round( np.log(0.5*Accuracy)/(np.log(RT_loss)) )
		
		# Calculate the field entering the cavity
		Ein = self.E_in.copy()
		Ein.change_n(new_n = self.I_input.n2)
		Ein = Ein.Ref_Trans_mirror(Iin = self.I_input)[1]
		
		Num_point_grid = self.E_in.Grid.Num_point
		self.Store_RT_field = np.zeros((Num_point_grid, Num_point_grid,num_iter), dtype = np.complex128) # where we will store the field after each round trip
		
		Field_circ = Ein.copy()
		
		ii = 0
		while ii <  num_iter:
			# print(self.Store_RT_field.shape)			
			# print(self.Store_RT_field.dtype)
			# print(Field_circ.shape)
			# print(Field_circ.dtype)
			
			self.Store_RT_field[:,:,ii] = Field_circ.Field
			Field_circ = self.do_round_trip(Field_circ)
			ii+=1
		
		# Now will sum the field according to the round trip phase
		Num_point_scan = 400
		self.Scan_RT_phase = np.linspace(0,2*np.pi,Num_point_scan)
		self.Scan_power = np.zeros_like(self.Scan_RT_phase)
		
		ii = 0
		while ii < Num_point_scan:
			Field_reconstructed = np.zeros((Num_point_grid, Num_point_grid), dtype = np.complex128)
			
			jj = 0
			while jj < num_iter:
				Field_reconstructed += self.Store_RT_field[:,:,jj] * np.exp(1j *jj * self.Scan_RT_phase[ii])
				jj += 1
			
			Field_circ.Field = Field_reconstructed
			self.Scan_power[ii] = Field_circ.power()
			ii += 1
		
		return self
		
		


class field(object):
	"""
	Class to represent the 2D electric field
	"""
	def __init__(self, grid, w = None, Rc = 1E99, w0 = None, z = 0, q = None, power = 1, mode = 'HG 0 0'):

		self.grid = grid
		
		self.ref_index = 1 # refractive index
		self.wavelength = 1064E-9 # [m]
        
        # TODO does this need to include ref_index??
		self.k_prop = 2 * np.pi / self.wavelength             
		self.mode = mode            
            
		if w is not None:
			q_start = 1/complex(1/Rc,-self.wavelength / (np.pi * w**2))
		elif w0 is not None:
			q_start = complex(z , np.pi * w0**2 / self.wavelength)
		else:
			q_start = q

		self.amplitude = np.exp(-1j*self.k_prop*self.grid.r_squared/(2*q_start))
		beam_radius = np.sqrt( 1/(-np.imag(1/q_start)*np.pi/(self.wavelength)) )

		mode, m, n = readModeName(self.mode)

		if mode == 'HG':
			# Prepate the Hermite polynomials                
			m_vec = np.zeros(m+1); m_vec[-1] = 1
			n_vec = np.zeros(n+1); n_vec[-1] = 1               

			self.amplitude *= np.polynomial.hermite.hermval(np.sqrt(2)/beam_radius * self.grid.X,m_vec)
			self.amplitude *= np.polynomial.hermite.hermval(np.sqrt(2)/beam_radius * self.grid.Y,n_vec)
       
		elif mode == 'LG':
			self.amplitude *= (2 * self.grid.r_squared / beam_radius**2) ** (np.abs(n)/2)
			self.amplitude *= sp.special.eval_genlaguerre(m, n, 2*self.grid.g_squared / beam_radius**2)
			self.amplitude *= np.exp(1j *n *np.arctan2(self.grid.Y,self.grid.X))

		else:
			print('field() problem: the mode name must be either HG or LG')
			return
            
		self.normalise(power = power)          

     
	def __add__(self,field2):  # overload the + operator
		out = self.copy()
		# Check if the two E_field could be added, it must have the same refractive index,
		# be defined on the same grid and have the same wavelength
		
		if self.grid != field2.grid:
			print('Could not add two fields defined on a different grid')
			return

		if self.ref_index != field2.ref_index:
			print('Could not add two fields with a different refractive index')
			return
			
		if self.wavelength != field2.wavelength:
			print('Could not add two fields with a different wavelength')
			return

		out.amplitude = self.amplitude + field2.amplitude	
		return out
	
	
	def __mul__(self,var_scalar):  # overload the * operator, mulitply the E_field with a scalar
		out = self.copy()

		out.amplitude = self.amplitude * var_scalar	
		return out

	def copy(self):
		out = copy.deepcopy(self)
		out.grid = self.grid
		return out
	
	def power(self):
		# Calculate the power in Watt of the field
		return np.sum(np.sum(np.abs(self.amplitude)**2)) * (self.grid.xstep * self.grid.ystep)
		
	def normalise(self, power = 1):
		if power == 0:
			self.amplitude = self.amplitude * 0
		else:
			current_power = self.power()
			if power != 0:
				self.amplitude = (self.amplitude / np.sqrt(current_power)) * np.sqrt(power)
		return self


	def plot(self):
		# Plot the 2D intensity of the field
		Z = np.abs(self.amplitude)**2

		fig = plt.figure()
		ax = plt.axes()
		                  
		xaxis_min = self.grid.xaxis.min()
		xaxis_max = self.grid.xaxis.max()
		yaxis_min = self.grid.yaxis.min()
		yaxis_max = self.grid.yaxis.max()
    
		ax.set_xlim(xaxis_min,xaxis_max)
		ax.set_ylim(yaxis_min,yaxis_max) 

		im = ax.imshow(Z,extent=[xaxis_min,xaxis_max,yaxis_min,yaxis_max])
		cb = fig.colorbar(im, ax=ax)

		
	def fitTEM00(self):
		# Fit the beam parameters for a fundamental Gaussian beam
		# For better accuracy, normalise the beam
		
		#self =  self.normalise();            

		# Do the fit only on the central part of the grid
		radius_fit = self.grid.length/4;

		# Extract only the central part of the grid
		grid_4_Fit = np.extract(self.grid.r < radius_fit,self.grid.r_squared)
		Power_4_Fit = np.extract(self.grid.r < radius_fit,np.abs(self.amplitude)**2)

		def Func_4_Fit(xdata, p): # The fitting function
			return  p[0]*(np.exp(-2*xdata/(p[1]**2)))         

		def Err_4_Fit(p, xdata, ydata): # The error function
			return  np.sum(abs(Func_4_Fit(xdata, p) - ydata)**2)
    
		res1 = minimize(Err_4_Fit, x0 = [np.max(Power_4_Fit), radius_fit], args=(grid_4_Fit, Power_4_Fit), method = 'Nelder-Mead')
		Radius_fitted = res1.x[1]            
		#print(Radius_fitted)            

		# First find an approximation for the RoC of the wavefront, do that in 1D
		Cross_sec_axis = np.extract(np.abs(self.grid.xaxis) < Radius_fitted*3,self.grid.xaxis)
		Cross_sec_phase_x = self.amplitude[self.grid.xpoints/2,:]
		Cross_sec_phase_x = np.extract(np.abs(self.grid.xaxis) < Radius_fitted*3,Cross_sec_phase_x)
		Cross_sec_phase_x = np.unwrap(np.angle(Cross_sec_phase_x))
    
		#plt.plot(Cross_sec_axis, Cross_sec_phase_x)

		pol_coeff = np.polyfit(Cross_sec_axis, Cross_sec_phase_x, 2)          
		if pol_coeff[0] != 0:
			beam_curvature_fit = - self.k_prop/(2*pol_coeff[0])
		else:
			beam_curvature_fit = 1E99
		
		#RoC_fitted = beam_curvature_fit
		#print(beam_curvature_fit)
		#options={'disp': True}               
 
		# # Do the fit on the complex curvature now
		grid_4_Fit = np.extract(self.grid.r < Radius_fitted*2,self.grid.r_squared)
		Amp_4_Fit = np.extract(self.grid.r < Radius_fitted*2,self.amplitude)

		def Err_4_Fit(p, xdata, ydata): # The error function
			Fct_return = p[0]*(np.exp(-2*xdata/(Radius_fitted**2) -1j*self.k_prop*xdata/(2*p[1]) + 1j*p[2]) )
			return  np.sum((np.abs(Fct_return - ydata))**2)
							
		res2 = minimize(Err_4_Fit, x0 = [np.max(abs(Amp_4_Fit)), beam_curvature_fit, 0.1], args=(grid_4_Fit, Amp_4_Fit), method = 'Nelder-Mead')
		RoC_fitted = res2.x[1]
		# #print(RoC_fitted)

		return ( Radius_fitted , RoC_fitted )


	def propagate(self , distance = 1):
		# Function to propagate the beam over a distance (simple FFT)
		
		plD = np.pi*self.wavelength*distance/self.ref_index
		propMatrix = np.exp(-1j*self.k_prop*distance) * np.exp(1j*plD*self.grid.fft_ir_squared)
		
		tmpfield = np.fft.fft2(self.amplitude)
		tmpfield = tmpfield * propMatrix
		self.amplitude = np.fft.ifft2(tmpfield)

		return self

        # TODO check that the above code returns the same as the old code below
	
		"""
		Propa_mat = np.exp( 1j*( -self.k_prop * dist + \
					np.pi*(self.Wavelength/self.Refractive_index) * \
					(self.grid.D2_FFT_X**2 + self.grid.D2_FFT_Y**2)* dist ))
           
		temp_Fourier_field = fftpack.fftshift(fftpack.fft2(self.amplitude))    
           
		self.amplitude = fftpack.ifft2(fftpack.ifftshift(temp_Fourier_field * Propa_mat ))
		return self
		"""
		
	def quickPropagate(self, propMatrix):
		# Useful if the propagation has already been pre-calculated
		tmpfield = np.fft.fft2(self.amplitude)
		tmpfield = tmpfield * propMatrix
		self.amplitude = np.fft.ifft2(tmpfield)
		return self

	def scalePropagate(self, distance, scale, newGrid):
		# FFT propagation code with an adaptive grid size.
		# Propagates to a scaled coordinate system, see Virgo Book of
		# Physics pages 179-184.
		# The scaling factor is given as w1 / w0 with w0 the beam size
		# at the start of propagation and w1 the expected beam size at
		# the end of propatation.
		#
		# NOT YET TESTED

		# For now, throw an error if distance == 0.
		# If distance is zero, the function should return a simple
		# scaled, unpropagated field.
		# Currently, if distance == 0, a division by zero occurs due to
		# the use of an inverse z0 function instead of simple z0 (which
		# is itself used to avoid scale == 1 causing a division by zero).
		assert(distance != 0)

		if scale <= 0:
			raise Exception('Specified scale factor must be > 0')

		if newGrid.xpoints != self.grid.xpoints or newGrid.ypoints != self.grid.ypoints:
			raise Exception('New grid must have same number of points as existing grid')

		plD = np.pi * self.wavelength * distance / (scale * self.ref_index)

		# invz0 is used to avoid division by zero if scale == 1
		invz0 = (scale - 1.0) / distance
	
		# initial scaling
		field = self.amplitude * np.exp(-1j * self.k_prop * self.grid.r_squared * invz0 / 2.0)
		field = np.fft.fft2(field)
		# scaled propagator
		field = field * np.exp(-1j * self.k_prop * distance) * np.exp(1j * plD * self.grid.fft_ir_squared)
		field = np.fft.ifft2(field)
		# final scaling
		self.amplitude = field * np.exp(1j * self.k_prop * newGrid.r_squared * invz0 * (1 + distance * invz0) / 2) / scale

		# update grid
		self.grid = newGrid

	def passAperture(self , diam = 0, shape = 'round'):
		# Pass through a round aperture of a fixed diameter
		# shape = 'round' or 'square'
		mask_aperture = np.zeros((self.grid.xpoints,self.grid.ypoints))            

		if shape == 'square':
			np.place(mask_aperture, np.logical_and(np.abs(self.grid.X) < diam, \
			abs(self.grid.Y) < diam), vals =  1)
     
		elif shape == 'round':
			np.place(mask_aperture, self.grid.r < diam/2, vals =  1)    
			
		self.amplitude = self.amplitude*mask_aperture
    
		return self           


	def calculateOverlap(field1,field2):
		# return the complex overlap integral between 2 instance of E_Field           
		E3 = field1.copy(); E4 = field2.copy()            
		E3.normalise()
		E4.normalise()
		overlap = np.sum(np.sum(E3.amplitude * np.conj(E4.amplitude))) * (field1.grid.xstep * field1.grid.ystep); 
		return overlap

		
	def expandHOM(Ein,Ebasis,max_mode_order):
		# Expand Ein on the basis defined by Ebasis (a TEM00)
		# The results is given as the normalised power of the Ein in the basis of the E_basis            
		result_vec = np.arange(max_mode_order+1)
		mat_overlap = np.zeros((max_mode_order+1,max_mode_order+1))

		# Fit the parameters of E_basis to define later the HOM
		rad_basis, curv_basis = Ebasis.fitTEM00()

		for ind_m in result_vec:               
			for ind_n in range(0,max_mode_order+1):
				mode_name = 'HG ' + str(ind_m) + ' ' + str(ind_n)
				Ebasis = field(Ein.grid, w = rad_basis, Rc = curv_basis, mode = mode_name)
				Coeff_overlap = field.calculateOverlap(Ein,Ebasis)                   
				mat_overlap[ind_m,ind_n] =  abs(Coeff_overlap)**2                   
				#print(mode_name,mat_overlap[ind_m,ind_n])
		
		result_vec2 = []
		for ind_m in np.arange(max_mode_order+1):
			tmp_sum = 0
			for ind_n in range(0,ind_m+1):
				tmp_sum += mat_overlap[ind_m-ind_n,ind_n]
			result_vec2.append(tmp_sum)
			
		#print(result_vec2)
		return result_vec2

	def change_n(self, new_n = None):	
		# Change the refractive index of the medium where the E_field is defined
		
		if new_n is None:
			print('change_n(): the new refractive index must be given')
			return
			
		if self.Refractive_index == new_n:
			print('change_n(): old and new refractive index are the same')
			return
			
		self.k_prop = self.k_prop * (new_n/self.ref_index)
		self.ref_index = new_n
		return self
		
		
	def transmission_lens(self, focal_length):
		# Transmission through a lens of a given focal length
		OPD_lens = ((2 * focal_length) - np.sign(focal_length)*np.sqrt((2*focal_length)**2 - self.grid.r_squared) )*2        
		WF_lens = np.exp(1j * self.k_prop * OPD_lens)

		self.amplitude = self.amplitude * WF_lens
		
			
	def ref_trans_mirror(self, Iin = None, RoC = None, Reflectivity = None):
		# Reflection and transmission from an interface or simple reflection from a mirror with a given curvature.
		# This function return 2 E_Fields, the reflected and the transmitted fields
		
		E_Ref = self.copy()
		E_Trans = self.copy()
		
		if Iin is not None: # Check first if an interface is given 
			# Check if the reflectivity has to be override (will not touch at the transmission) 
			if Reflectivity is not None:
				ref_interface = np.sqrt(Reflectivity)
			else:
				ref_interface = Iin.r
		
			if self.Refractive_index == Iin.n1:
				# Calculation of the wavefront distortion
				WF_Mirror_trans = np.exp(1j * E_Trans.k_prop * ((Iin.n2 - Iin.n1)/ Iin.n1) * Iin.surface) * Iin.mask * Iin.t
				E_Trans.change_n(Iin.n2)

				WF_Mirror_ref = np.exp(-1j * E_Ref.k_prop * 2 * Iin.surface) * Iin.mask * ref_interface
			
			elif self.Refractive_index == Iin.n2:
				# Calculation of the wavefront distortion
				WF_Mirror_trans = np.exp(- 1j * E_Trans.k_prop * ((Iin.n1 - Iin.n2)/ Iin.n2) * np.fliplr(Iin.surface)) * \
				np.fliplr(Iin.mask) * Iin.t
				E_Trans.change_n(Iin.n1)
				
				WF_Mirror_ref = np.exp(1j * E_Ref.k_prop * 2 * np.fliplr(Iin.surface)) * np.fliplr(Iin.mask) * ref_interface
					
			else:
				print('Ref_Trans_mirror(): refractive index not matching')
		
		elif RoC is not None: # The user give a radius of curvature for the mirror
			if Reflectivity is None:
				ref_interface = 1
			else:
				ref_interface = np.sqrt(Reflectivity)
			
			# Calculation of the wavefront distortion in reflection
			Surface_Mirror = (RofC - np.sign(RofC) * np.sqrt(RofC**2 - self.grid.r_squared))*2
			WF_Mirror_ref = np.exp(1j * E_Ref.k_prop * Surface_Mirror) * ref_interface
    		
			# Do a dummy one for the transmission
			WF_Mirror_trans = np.ones((self.grid.Num_point,self.grid.Num_point))
			
		else:
			print('Ref_Trans_mirror(): error with the input arguments')
			return
		
		# Applied the wavefront distortion to the field	
		E_Trans.amplitude = E_Trans.amplitude * WF_Mirror_trans
                    
		E_Ref.amplitude = E_Ref.amplitude * WF_Mirror_ref
		E_Ref.amplitude = np.fliplr(E_Ref.amplitude)
		
		return E_Ref,E_Trans
		
	def quick_reflection(self, WF_ref): # Useful if the wavefront distortion has already been pre-calculated	
		self.amplitude *= WF_ref
		self.amplitude = np.fliplr(self.amplitude)
		
		
		
		


class interface(object):
	"""
	Class to represent the 2D wavefront distortion (or surface)
	"""
	
	def __init__(self, _grid, _RoC = 1E99, _CA = 1E99, _T = 0.1, _L = 0, _n1 = 1, _n2 = 1.45):
		
		self.grid = _grid
		
		if _RoC == 0:
			_RoC = 1E99

		self.surface =  -(_RoC - np.sign(_RoC) * np.sqrt(_RoC**2 - self.grid.r_squared) )
	
		self.n1 = _n1
		self.n2 = _n2
	
		self.T = _T
		self.L = _L
		
		self.t = 1j*np.sqrt(self.T)
		self.r = np.sqrt(1 - (self.T + self.L))
		
		# Create the mask for the aperture
		self.mask = np.zeros((self.grid.xpoints,self.grid.ypoints))
		np.place(self.mask, self.grid.r < _CA/2, vals =  1)		
	
	
	def plot(self):
		"""
		Plotting the wavefront distortion
		"""
		fig, ax = plt.subplots()
                  
		xaxis_min = self.grid.xaxis.min()
		xaxis_max = self.grid.xaxis.max()
		yaxis_min = self.grid.xaxis.min()
		yaxis_max = self.grid.xaxis.max()
    
		ax.set_xlim(xaxis_min,xaxis_max)
		ax.set_ylim(yaxis_min,yaxis_max)     

		im = ax.imshow(self.surface,extent=[xaxis_min,xaxis_max,yaxis_min,yaxis_max])
		cb = fig.colorbar(im, ax=ax)

class grid():
	"""
	Data structure to describe the size and axes for a (x,y) data array
	of complex beam amplitudes. Also contain also data structures for
	FFT propagation
	"""
	def __init__ (self, xpoints, ypoints, xlength, ylength, xoffset=0, yoffset=0):

		self.xpoints=xpoints # [number of tiles]
		self.ypoints=ypoints # [number of tiles]
		self.xlength=xlength # [m]
		self.ylength=ylength # [m]
		self.xoffset=xoffset # [m]
		self.yoffset=yoffset # [m]

		self.length = (xlength + ylength )/2.0 #[m]
		
		# compute x and y axis
		self.xstep=self.xlength/self.xpoints # [m]
		self.ystep=self.ylength/self.ypoints # [m]
		xvector= np.arange(self.xpoints)
		yvector= np.arange(self.ypoints)
		self.xaxis=-self.xlength/2.0 + self.xstep/2.0 + xvector*self.xstep + self.xoffset
		self.yaxis=-self.ylength/2.0 + self.ystep/2.0 + yvector*self.ystep + self.yoffset

		# and some useful variables based on the axis
		self.X,self.Y = np.meshgrid(self.xaxis,self.yaxis)
		self.r_squared = (self.X)**2 + (self.Y)**2
		self.r = np.sqrt(self.r_squared)
		self.angle = np.arctan2(self.Y,self.X)

		# compute frequency axis
		self.xaxis_fft = np.fft.fftshift(np.fft.fftfreq(self.xpoints))/self.xstep
		self.yaxis_fft = np.fft.fftshift(np.fft.fftfreq(self.ypoints))/self.ystep

		# some useful variables based on the frequency axis
		self.fft_X,self.fft_Y = np.meshgrid(self.xaxis_fft, self.yaxis_fft)
		self.fft_ir_squared= np.fft.ifftshift((self.fft_X)**2+(self.fft_Y)**2)


def readModeName(_str_name):
	"""
	Utility function: Take the full mode name 'LG 2 3' and return the
	mode mode name 'LG' and the mode numbers (2,3) 
	"""
	family = _str_name[0:2]
	mode_number = _str_name[3:] 
	
	# Slice the mode number
	m, n = mode_number.split()
	
	return family, int(m), int(n)
	

	
