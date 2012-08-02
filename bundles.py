'''Spectroscopy Bundles Submodule
	Author: Alex M. Pronschinske
	
	List of classes:
		SpecBundle
		IVSTSBundle
		ZVSTSBundle
	List of functions: -none-
'''

# built-in modules
import os
import re
import time
import warnings

# third-party modules
import numpy as np
from random import sample as randsample
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# internal modules
import spectroscopy.numerical as specNum
import spectroscopy.specPlot as specPlot
import spectroscopy.specFit as specFit
import spectroscopy.specTools as specTools
import spectroscopy.importers as imps

#===============================================================================
class SpecBundle(object):
	'''Spectra Bundle Class
	
	This base class will enable specific processing methods for a list of
	ndarrays representing data curves.  Since SpecBundle is a base class it
	cannot be instantiated, only sub-classes can be instantiated.  An instance
	can be instantiated from data in memory or from a data file.
	
	Instantiation Args:
		* There are multiple input patterns, see examples.
		Acceptable keyword arguments:
			file_name (str)
			file_type (str)
	Instance Attributes:
		file_name (str): Data file of origin
		X (ndarray): Common X array
		allY (list): List of ndarray's containing individual np.ctra
		N (int): Number of np.ctra in the bundle
		curve_len (int): Number of points in any given np.ctra
		units (list): Length 2 list of str's containing unit labels for X & Y
		has_ramp_reversal (bool)
	Class Methods:
		__add__
		__getitem__
		__iadd__
		__isub__
		__iter__
		__setitem__
		__sub__
		add_spec
		extend
		copy
		insert_spec
		pop
		reverse
		coavg
		correl_matrix
		deriv_sg
		fan
		find_peaks
		plot_density
		ramp_reversal_detect
		ramp_reversal_fix
		save_multispec
	Examples:
		allspec = ZVSTSBundle('data.zvms.asc', file_type='matrix_asc')
		allspec = ZVSTSBundle(V, allZ, ('V', 'nm'), file_name='data.zvms.asc')
		allspec = IVSTSBundle(V, allI, ('V', 'pA'), file_name='data.ivms.asc')
	'''
	
	def __init__(self, *args, **kwargs):
		if type(self) is SpecBundle:
			raise RuntimeError(
				'SpecBundle is a base class and cannot be instatiated'
			)
		# END if
		
		# For one string argument assumed it is a file path
		# Call _import_file method of the sub-class for reading the file in
		if len(args) == 1 and type(args[0]) is str:
			self.file_name = args[0]
			if 'file_type' in kwargs:
				self._import_file(self.file_name, kwargs['file_type'])
			else:
				self._import_file(self.file_name)
			# END if
		else:
			# TODO: add type validation
			self.X = args[0]
			self.allY = args[1]
			if len(args) >= 3 and type(args[2]) is str:
				self.units = args[2]
			else:
				self.units = None
			# END if
			
			if 'file_name' in kwargs.keys():
				self.file_name = kwargs['file_name']
			else:
				self.file_name = None
				#raise Exception('SpecBundle must have file_name arg defined')
			# END if
		# END if
		
		self.__rr = None
	# END __init__
	
	# Property descriptor functions
	#------------------------------
	@property
	def N(self):
		return len(self.allY)
	# END N
	
	@property
	def curve_len(self):
		return len(self.X)
	# END curve_len
	
	@property
	def has_ramp_reversal(self):
		if self.__rr is None: self.__rr = self.ramp_reversal_detect()
		return self.__rr
	# END has_ramp_reversal
	
	# Magic methods
	#--------------
	
	def __add__(self, other_bun):
		bun_sum = self.copy
		bun_sum.extend(other_bun)
		return bun_sum
	# END __add__
	
	def __getitem__(self, key):
		return self.allY[key]
	# END __getitem__
	
	def __len__(self):
		return len(self.allY)
	# END __len__
	
	def __iadd__(self, other_bun):
		# TODO: add polymorphism so that a single spectrum or list of spectra
		#       can be added as well
		self.extend(other_bun)
	# END __iadd__
	
	def __isub__(self, arg):
		if type(arg) is int:
			self.pop(arg)
		elif type(arg) is list:
			for i in list: self.pop(i)
		# END if
	# END __isub__
	
	def __iter__(self):
		'''SpecBundle will iterate over the allY attribue'''
		for i in range(self.N): yield self[i]
	# END __iter__
	
	def __setitem__(self, key, value):
		if type(value) == list:
			value = np.array(value)
		elif type(value).__name__ == 'ndarray':
			pass
		else:
			raise TypeError(
				'Items in SpecBundle must be of type list or ndarray not '
				+ type(value).__name__
			)
		# END if
		self.allY[key] = value
	# END __setitem__
	
	def __sub__(self, arg):
		new_bun = self.copy()
		if type(arg) is int:
			new_bun.pop(arg)
		elif type(arg) is list:
			for i in list: new_bun.pop(i)
		# END if
		return new_bun
	# END __sub__
	
	# List-like methods
	#------------------
	def add_spec(self, *args):
		'''Adds a single spectra (ndarray) to the end of self.allY
		
		The added spectrum is assumed to be the same length.
		
		Examples:
			bun.add_spec(Y_new)
		'''
		
		for Y in args: self.allY.append(Y)
	# END add_spec
	
	def extend(self, *args):
		'''Adds all spectra in a SpecBundle to self
		
		The added spectra are assumed to be the same length.  This method is
		synonymous with the "+=" operator.
		
		Examples:
			bun_1.extend(bun_2, bun_3, bun_4)
		'''
		
		for bundle in args:
			self.allY.extend(bundle.allY)
			self.file_name += ',{}'.format(bundle.file_name)
	# END extend
	
	def copy(self):
		'''Make a copy of self'''
		
		return type(self)(
			self.X, self.allY, self.units, file_name=self.file_name
		)
	# END copy
	
	def insert_spec(self, index, Y):
		'''Adds a single spectra (ndarray) anywhere into self.allY
		
		The added spectrum is assumed to be the same length.
		
		Examples:
			bun.insert_spec(3, Y_new)
		'''
		
		self.allY.insert(index, Y)
	# END insert_spec
	
	def pop(self, index=0):
		'''Removes and returns a spectrum from the specified index
		
		Args:
			index = 0 (int)
		Examples:
			popped_Y = bun.pop(3)
		'''
		return self.allY.pop(index)
	# END pop
	
	def reverse(self):
		'''Reverses the self.allY list'''
		
		self.allY = self.allY[::-1]
	# END reverse
	
	# SpecBundle-specific methods
	#----------------------------
	def coavg(self, output_stdev=False):
		'''Co-average spectra
		
		The returned spectrum will be an average of all spectra in the bundle.
		
		Args:
			output_stdev = False (bool): If true the tuple, (mY, eY) is output
										instead of just the ndarray mY.
		Returns:
			(ndarray|tuple) Outputs the ndarray mY ("mean Y") and also
			optionally eY ("error in Y") which is the standard deviation.
		'''
		
		mY = np.zeros(self.curve_len)
		eY = np.zeros(self.curve_len)
		for Y in self.allY:
			mY = mY + Y
			eY = eY + Y**2
		# END for
		mY = mY / self.N
		eY = np.sqrt( eY/self.N - mY**2 )
		
		if output_stdev: return mY, eY
		else: return mY
	# END coavg
	
	def correl_matrix(self):
		'''Method for Finding the Correlation Matrix of the Bundle
		
		Args: none
		Returns:
			(matrix) A numpy matrix object that contains the correlations
			between all of the spectra in the bundle
		'''
		
		cmat = np.matrix(
			[ [1.0 for j in range(self.N)] for i in range(self.N) ]
		)
		
		for i in range(self.N):
			Yi = self.allY[i]
			for j in range(i+1, self.N):
				r = pearsonr( Yi, self.allY[j] )[0]
				r = np.log(r) - np.log(1.0-r)
				cmat[i,j] = r
				cmat[j,i] = r
			# END for
		# END for
		
		return cmat
	# END correl_matrix
	
	def deriv_sg(self, x_window, poly_order, deriv_order=0):
		'''Differential all spectra in bundle
		
		Calculates the derivative of the individual spectra using a Savitzky-
		Golay filter.  This is NOT an in-place method.
		
		Args:
			x_window (float): Width of the fitting window is units of the x-axis
			poly_order (int): Order of the fitting polynomial
			deriv_order (int): Order of the derivative
		Returns:
			(same type as self) SpecBundle sub-class where all spectra are
			corresponding derivatives of this bundle's spectra
		Examples:
			deriv_bun = bun.deriv_sg(0.5, 3, 1)
		'''
		
		# Validate input
		try:
			poly_order = np.abs(int(poly_order))
		except ValueError:
			raise ValueError("polynomial order has to be of int type")
		# END try
		
		dx = self.X[1] - self.X[0]
		i_window = int(x_window/dx)
		i_window += 1 - i_window%2
		if i_window < poly_order + 2:
			raise TypeError(
				"window is too small for the polynomials order"
				)
		# END if
	
		order_range = range(poly_order+1)
		half_window = (i_window -1) / 2
	
		# precompute coefficients
		b = np.mat(
			[
				[k**i for i in order_range] for k in range(
					-half_window, half_window+1
				)
			]
		)
		m = np.linalg.pinv(b).A[deriv_order]
		
		dYs = []
		for Y in self.allY:
			# pad the function at the ends with reflections
			left_pad = Y[0] - np.abs( Y[1:half_window+1][::-1] - Y[0] )
			right_pad = Y[-1] + np.abs(Y[-half_window-1:-1][::-1] - Y[-1])
			Y = np.concatenate((left_pad, Y, right_pad))
			
			dY = np.power(-1, deriv_order) * np.convolve( m, Y, mode='valid')
			dY /= np.power(dx, deriv_order)
			dYs.append(dY)
		# END for
		
		return type(self)(self.X, dYs, file_name=self.file_name)
	# END deriv_sg
	
	def fan(self, delta, in_place=False):
		'''Spread out the spectra vertically
		
		Adds amount i*delta to each Y = self.allY[i].  This method is intended
		to be used on a bundle before plotting.
		
		Args:
			delta (float): Space between adjacent spectra in units of y-axis
			in_place = False (bool): Controls if the fanning is done in-place
		Returns:
			(same as self) Bundle of fanned-out spectra
		
		TODO: make this an in-place ONLY version of a function in the plotting
		      sub-module
		'''
		if in_place:
			for i in range(self.N):
				self[i] = self[i] + i*delta
			# END for
			return None
		else:
			out_specs = self.allY
			for i in range(self.N):
				out_specs[i] = out_specs[i] + i*delta
			# END for
			return type(self)(
				self.X, out_specs, self.units, file_name=self.file_name
			)
		# END if
	# END fan
	
	def find_peaks(self, xrng, A_min=0.0, full_output=False, write_all=False):
		'''Peak location analysis method
		
			Fits a gaussian peak to each spectrum in the bundle and returns
			information about the fit.  The gaussian function has the following
			form, y(x) = A * exp( -ln(16) * ((x-x0)/w)**2 ) + y0
		
			Args:
				xrng (list): Peak search domain. "x-range"= [x_min, x_max]
				A_min (float): Minimum acceptable peak amplitude
				full_output (bool): Will populate return array with dict's
					describing each peak.
				write_all (bool): Will write a detailed plot of each fit to the
					scratch directory
			Returns:
				(list) List of peak locations
			Returns w/ full_output=True:
				(list) List of dictionaries with key values:
					'i': spectrum index in bundle,
					'x0': peak location,
					'w': peak full-width at half-maximum,
					'A': peak amplitude (ie. A=ymax-y0),
					'y0': peak baseline (ie. y value at +/- inf),
					'Rsq': R**2 value for the goodness of the fit
		'''
		
		if xrng[1] < xrng[0]: xrng.reverse()
		
		peaks = []
		accepted_Y = []
		accepted_i = []
		allX0 = []
		allW = []
		ierr = 0
		#errMax = int(0.55*len(self.allY))
		
		# translate x range to an index range
		irng = [None, None]
		gamma = 1.5
		a = np.mean(xrng) - gamma*(xrng[1]-xrng[0])/2.0
		b = np.mean(xrng) + gamma*(xrng[1]-xrng[0])/2.0
		# TODO: this for loop is not ideal, find i's faster by calculating them
		for i, x in enumerate(self.X):
			if a <= x and x < b and irng[0] is None:
				irng[0] = i
			elif b <= x and irng[1] is None:
				irng[1] = i
		# END for
		if irng[1] is None: irng[1] = len(self.X)-1
		
		# Create a folder in the scratch directory for "write_all" files
		if write_all:
			saveDir = (
				'/home/alex/scratch/'
				+ os.path.basename(self.file_name)[:-4] + '/'
			)
			if not os.path.exists(saveDir):
				os.mkdir(saveDir)
			else:
				for f in os.listdir(saveDir): os.remove(saveDir+f)
			# END if
		# END if
		
		# Working loop
		for i, Y in enumerate(self.allY):
			rejectMsg = None
			
			# Fit the spectrum
			[x0, sig, A, y0, Rsq] = specFit.gaussFit(
				self.X[irng[0]:irng[1]], Y[irng[0]:irng[1]]
			)
			w = np.sqrt(8.0*np.log(2.0))*sig
			
			if A <= 0:
				ierr += 1
				rejectMsg = 'Negative Peak'
			elif A < A_min:
				ierr += 1
				rejectMsg = 'Peak too Small: A={0:0.3f}'.format(A)
			elif A/w > 18.0:
				ierr += 1
				rejectMsg = 'Peak aspect ratio too high'
			elif x0 < xrng[0] or xrng[1] < x0:
				ierr += 1
				rejectMsg = 'Peak Outside X Bounds: x0={0:0.2f}'.format(x0)
			#elif (1-Rsq) > 2.0E-3:
			#	ierr += 1
			#	rejectMsg = 'Low R^2: (1-R^2)= {0:0.2e}'.format(1-Rsq)
			elif 2.0 < w:
				ierr += 1
				rejectMsg = 'Peak Outside W Bounds: w={0:0.2f}'.format(w)
			else:
				peaks.append(
					{'i': i, 'x0': x0, 'w': w, 'A': A, 'y0': y0, 'Rsq': Rsq}
				)
				accepted_Y.append(Y)
				accepted_i.append(i)
				allX0.append(x0)
				allW.append(w)
			
			if write_all:
				Xf = self.X[irng[0]:irng[1]]
				Yf = A * np.exp( -0.5*((Xf-x0)/sig)**2 ) + y0
				Xextra = np.concatenate((self.X[:irng[0]], self.X[irng[1]:]))
				Yextra = A * np.exp( -0.5*((Xextra-x0)/sig)**2 ) + y0
				if w < 4.0 and A < 12:
					plt.plot(
						self.X, Y, '.',
						Xf, Y[irng[0]:irng[1]],
						Xf, Yf, Xextra, Yextra, '.',
						[x0, x0], [y0, y0+A],
						[x0-0.5*w, x0+0.5*w], [y0+0.5*A, y0+0.5*A]
					)
				else:
					plt.plot(
						self.X, Y, '.',
						Xf, Y[irng[0]:irng[1]],
						Xf, Yf
					)
				plt.grid(True)
				titleStr = (
					'i= {0:02d}, x0={1:0.3f}, w={2:0.3f}, A={3:0.2f}, ' +
					'1-R^2= {4:0.2e}'
				)
				titleStr = titleStr.format(i, x0, w, A, 1-Rsq)
				if rejectMsg is None:
					plt.title(titleStr)
					plt.savefig(saveDir + 'goodfit_{0}.png'.format(len(allX0)))
				else:
					titleStr = rejectMsg + '\n' + titleStr
					plt.title(titleStr)
					plt.savefig(saveDir + 'badfit_{0}.png'.format(ierr))
				plt.close()
			# END if
			
			#if ierr == errMax and write_all:
			#	print 'error limit of {0}/{1} reached!'.format(errMax, self.N)
		#END for
		
		if full_output:
			return peaks
		else:
			return allX0
		# END if
	# END find_peaks
	
	def plot_density(self, ax=None, mkrcolor='blue'):
		'''Calls specPlot.plot_density with self'''
		return specPlot.plot_density(self, ax, mkrcolor)
	# END plot_density
	
	def ramp_reversal_detect(self, N=3):
		'''Detect horizontally concatenated (i.e. ramp-reversed) data
			
			This method will sample N random spectra from the bundle and look
			for a large jump in the data at the halfway point as a marker of
			concatenated forward and backward spectrum sweeps.
			
			Args:
				N = 3 (int): Number of spectra to sample for detection
			Returns:
				(bool)
		'''
		
		if N > self.N:
			someY = self.allY
		else:
			someY = randsample(self.allY, 3)
		# END if
		
		dratio = 0
		L = self.curve_len
		if L%2 == 0:
			for Y in someY:
				dmid = Y[L/2-1] - Y[L/2]
				dwhole = Y[-1] - Y[0]
				dratio += abs(dmid/dwhole)/N
		else:
			dratio = 0
		# END if
		
		if dratio > 0.5:
			return True
		else:
			return False
		# END for
	# END ramp_reversal_detect
	
	def ramp_reversal_fix(self, opt='double'):
		'''Separate horizontally concatenated (i.e. ramp-reversed) data
			
			Args:
				opt (str): Name of method for handling separated data. Options
					include: double (one bundle w/ both forward & backward),
					average (one bundle w/ forward-backward averages), and
					split (two bundles, a forward bundle & a backward bundle)
			Returns:
				(SpecBundle) or (SpecBundle, SpecBundle)
		'''
		
		newX = np.linspace(self.X[0], self.X[-1], self.curve_len/2)
		
		if opt == 'double':
			doubY = []
			h = self.curve_len/2
			for Y in self.allY:
				doubY.append(Y[:h])
				doubY.append(Y[h:])
			# END for
			return type(self)(newX, doubY, self.units, file_name=self.file_name)
		elif opt == 'average':
			avgY = []
			h = self.curve_len/2
			for Y in self.allY:
				avgY.append( (Y[:h]+Y[h:])/2 )
			# END for
			return type(self)(newX, avgY, self.units, file_name=self.file_name)
		elif opt == 'split':
			fY = []
			rY = []
			h = self.curve_len/2
			for Y in self.allY:
				fY.append(Y[:h])
				rY.append(Y[h:])
			# END for
			forward_bndl = type(self)(
				newX, fY, self.units, file_name=self.file_name
			)
			reverse_bndl = type(self)(
				newX, rY, self.units, file_name=self.file_name
			)
			return forward_bndl, reverse_bndl
		else:
			raise ValueError(
				'valid opt values include: "double", "average", or "split"'
			)
		# END if
	# END ramp_reversal_fix
	
	def save_bundle(self, file_path):
		'''Save self as a text file
		
		This method will write the data to a text file with the following
		format:
		# <header lines>
		
		<col units>   <col units>   ...
		<col axis>    <col axis>    ...
		<data point>  <data point>  ...
		
		Args:
			file_path (str): save location and file name with appropriate
							 file extension
		'''
		
		f = open(file_path, 'w')
		
		# Write header
		f.write('# File Format = ASCII\n')
		f.write('# Source files:\n')
		files = re.split(r',', self.file_name)
		for file_name in files:
			file_name = re.sub(r'^(.*?)([^/]+)$', r'\2', file_name)
			f.write('#  {}\n'.format(file_name))
		f.write(
			'# Created: ' + time.strftime(r'%m/%d/%Y %H:%M:%S') + '\n'
		)
		f.write('\n')
		
		# Write column labels
		units_ln = self.N*'{0:12}\t{1:12}\t'.format(*self.units) + '\n'
		f.write(units_ln)
		axes_ln = self.N*'X           \tY           \t' + '\n'
		f.write(axes_ln)
		
		# Write data
		for i_ln in range(self.curve_len):
			data_ln = ''
			for Y in self:
				data_ln += '{0:0<+8f}\t{1:0<+8f}\t'.format(
					self.X[i_ln], Y[i_ln]
				)
			# END for
			data_ln += '\n'
			f.write(data_ln)
		# END for
		
		f.close()
	# END save_bundle
# END SpecBundle

#===============================================================================
class IVSTSBundle(SpecBundle):
	'''I(V)-STS Spectra Bundle Class
	
	This SpecBundle sub-class is specifically for analyzing I(V)-STS data.
	
	Class Methods:
		_import_file
		zero_y0
		norm_deriv
		filter_nd_abs
	'''
	
	def _import_file(self, file_name, file_type=''):
		'''Import file method
		
		This method will assume the data is MATRIX data saved as an ASCII text
		file by SPIP and call that importer module, unless file_type is
		otherwise specified.
		'''
		
		if file_type is '':
			# Detect the file type
			file_type = 'matrix_asc'
		# END if
		
		if file_type is 'matrix_asc':
			X, allY, units = imps.import_matrix_asc_iv(file_name)
		elif file_type is 'scala':
			X, allY, units = imps.import_scala_iv(file_name)
		else:
			raise RuntimeError(
				'"{}" is an unknown file type'.format(file_name)
			)
		# END if
		
		# **Set object properties**
		self.X = X
		self.allY = allY
		self.units = units
		self.zero_y0()
	# END _import_file
	
	def zero_y0(self):
		'''Shifts all spectra vertically so they intersect (0, 0)
		'''
		
		X = self.X
		dx = X[1] - X[0]
		i0_est = int(-X[0]/dx)
		if X[i0_est] != 0 and X[i0_est+1] != 0:
			i0 = i0_est
			for j, Y in enumerate(self.allY):
				dy = Y[i0+1]-Y[i0]
				y0 = Y[i0] - (dy/dx)*X[i0]
				self.allY[j] = Y - y0
			# END for
		else:
			i0 = i0_est
			if X[i0_est+1] == 0: i0 += 1
			for j, Y in enumerate(self.allY):
				self.allY[j] = Y - Y[i0]
		# END if
	# END zero_y0
	
	def norm_deriv(self, x_window, poly_order):
		''' Calculate the normalized dI/dV (i.e. (V/I)*dI/dV) for all spectra
			
			This method does NOT operate in-place.
			
			Args:
				x_window (float): Differentiation parameter passed to deriv_sg
				plot_order (int): (see above)
			Returns:
				(SpecBundle) All spectra are corresponding normalized derivative
				of self's spectra
			Examples:
				ndivsts = ivsts.norm_deriv(0.1, 1)
		'''
		
		dYs = self.deriv_sg(x_window, poly_order, 1)
		sYs = self.deriv_sg(x_window, poly_order, 0)
		sYs.zero_y0()
		ndYs = dYs.copy()
		for i in range(ndYs.N):
			ndYs[i] = self.X * dYs[i] / sYs[i]
		
		return ndYs
	# END norm_deriv
	
	def filter_nd_abs(
		self, flvl, x_window=0.1, poly_order=1, std_out=False
	):
		'''Filters out data with missing points and relatively large spikes
		
		This method behaves like the pop method, where the popped elements are
		rejected spectra and self is left with only acceptable spectra.
		
		Args:
			flvl (float): Absolute rejection level in units of y-axis
			x_window = 0.1 (float)
			poly_order = 1 (int)
			std_out = False (bool): Controls print statements
		Returns:
			(same as self) Bundle containing the rejected spectra
		Examples:
			rejects = bun.filter_nd_abs(5)
			rejects = bun.filter_nd_abs(5, 0.5, 3, True)
		'''
		
		warnings.simplefilter('error', RuntimeWarning)
		
		X = self.X
		N = self.N
		rejects = []
		n = 0
		while n < N:
			Y = self.pop()
			# Assume the spectra is ramp-reversed data
			X_, Y_f, Y_r = specTools.split_spectrum(X, Y)
			
			missing_point = False
			try:
				ndY_f = specNum.norm_deriv(X_, Y_f, x_window, poly_order)
				ndY_r = specNum.norm_deriv(X_, Y_r, x_window, poly_order)
				m = abs(max( [ndY_f.max(), ndY_r.max()] ))
			except RuntimeWarning:
				m = 0
				missing_point = True
			# END try
			
			if m > flvl or missing_point:
				if std_out: print 'curve {0:03d}, m = {1}'.format(n, m)
				rejects.append(Y)
			else:
				self.add_spec(Y)
			# END if
			
			n += 1
		# END while
		
		return type(self)(self.X, rejects, self.units, file_name=self.file_name)
	# END filter_nd_abs
# END IVSTSBundle

#===============================================================================
class ZVSTSBundle(SpecBundle):
	'''z(V)-STS Spectra Bundle Class
	
	This SpecBundle sub-class is specifically for analyzing z(V)-STS data.
	
	Class Methods:
		_import_file
		make_relative
	'''
	
	def _import_file(self, file_name, file_type=''):
		'''Import file method
		
		This method will assume the data is MATRIX data saved as an ASCII text
		file by SPIP and call that importer module, unless file_type is
		otherwise specified.
		'''
		
		if file_type is '':
			# Detect the file type
			file_type = 'matrix_asc'
		# END if
		
		if file_type is 'matrix_asc':
			X, allY, units = imps.import_matrix_asc_zv(file_name)
		else:
			raise RuntimeError(
				'"{}" is an unknown file type'.format(file_name)
			)
		# END if
		
		# **Set object properties**
		self.X = X
		self.allY = allY
		self.units = units
		self.make_relative()
	# END _import_file
	
	def make_relative(self):
		'''Set all spectra to start from z=0'''
		X = self.X
		if X[0] < 0 and 0 < X[-1]:
			# find minimum and set it to zero
			raise RuntimeError('cannont handle z(V)-STS through V=0, yet...')
		elif 0 < X[0]:
			i0 = X[0]
		elif X[-1] < 0:
			i0 = X[-1]
		else:
			raise RuntimeError(
				'X list from "{}" may be backwards, dx={1:0.2f}'.format(
					self.file_name, X[1]-X[0]
				)
			)
		# END if
		
		for j, Y in enumerate(self.allY):
			self.allY[j] = Y - Y[i0]
		# END if
	# END make_relative
# END ZVSTSBundle