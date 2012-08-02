'''Spectroscopy Package Data Fitting Module
	Author: Alex M. Pronschinske
	
	List of classes: -none-
	List of functions:
		gaussFit
		gaussLinFit
		lorentzFit
		linfit
'''

# third-party modules
import numpy as np
from scipy.optimize import leastsq

#===============================================================================
def gaussFit(X, Y):
	'''Gaussian least-squares fit
	
	This function fits a gaussian function to the given data, it is assumed that
	the given data only contains ONE peak.  The functional form of the gaussian
	peak is as follows:
		y(x) = A * exp( -0.5 * ((x-x0)/sig)**2 ) + y0,
	where A is the height of the peak from the baseline (y0), x0 is the peak
	center, and sig is the standard deviation.  The initial parameter guess are:
		x0, y-weighted mean of the x values,
		sig, y-weighted standard deviation of the x values,
		A, maximum y value,
		y0, median value of the y data.
	This function also returns the R**2 "goodness" of the fit (Rsqr).
	
	Args:
		X (list-like): x-axis data
		Y (list-like): y-axis data
	Returns:
		(numpy.ndarray) [x0, sig, A, y0, Rsqr]
	'''
	
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	pInit = np.array([0.0, 0.0, max(Y), np.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = np.sqrt( ex1 - pInit[0]**2 )
	
	def GaussResid(params):
		x0 = params[0]
		sig = params[1]
		A = params[2]
		y0 = params[3]
		F = A * np.exp(-0.5*((X-x0)/sig)**2 ) + y0
		return F-Y
	
	(fitParams, cov_x, infodict, mesg, intFlag) = leastsq(
		GaussResid, pInit, full_output=True
		)
	if fitParams[1] < 0:
		fitParams[1] = abs(fitParams[1])
	
	Yf = infodict['fvec']
	Rsq = 1.0 - ( np.std(Y)**2/sum((Y-Yf)**2) )
	
	return np.append(fitParams, Rsq)
# END gaussFit

#===============================================================================
def gaussLinFit(X, Y):
	'''Gaussian plus linear background least-squares fit
	
	This function fits a gaussian function with a linear background to the given
	data, it is assumed that the given data only contains ONE peak.  The
	functional form is as follows:
		y(x) = A * exp( -0.5 * ((x-x0)/sig)**2 ) + y1*x + y0,
	where A is the height of the peak from the baseline (y1*x + y0), x0 is the
	peak center, and sig is the standard deviation.  The initial parameter
	guess are:
		x0, y-weighted mean of the x values,
		sig, y-weighted standard deviation of the x values,
		A, maximum y value,
		y1, zero
		y0, median value of the y data.
	This function also returns the R**2 "goodness" of the fit (Rsqr).
	
	Args:
		X (list-like): x-axis data
		Y (list-like): y-axis data
	Returns:
		(numpy.ndarray) [x0, sig, A, y1, y0, Rsqr]
	'''
	
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	#Parameter Array: [x0, sig, A, y1, y0]
	pInit = np.array([0.0, 0.0, max(Y), 0.0, np.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = np.sqrt( ex1 - pInit[0]**2 )
	
	def GaussResid(params):
		x0 = params[0]
		sig = params[1]
		A = params[2]
		y1 = params[3]
		y0 = params[4]
		F = A * np.exp(-0.5*((X-x0)/sig)**2 ) + y1*X + y0
		return F-Y
	
	(fitParams, cov_x, infodict, mesg, intFlag) = leastsq(
		GaussResid, pInit, full_output=True
		)
	if fitParams[1] < 0:
		fitParams[1] = abs(fitParams[1])
	
	Yf = infodict['fvec']
	Rsq = 1.0 - ( np.std(Y)**2/sum((Y-Yf)**2) )
	
	return np.append(fitParams, Rsq)
# END gaussFit

#===============================================================================
def lorentzFit(X, Y):
	'''Lorentzian least-squares fit
	
	This function fits a Lorentzian function to the given data, it is assumed
	that the given data only contains ONE peak.  The functional form of the
	Lorentzian peak is as follows:
		y(x) = A / (1.0 + 4.0*((x-x0)/w)**2 ) + y0,
	where A is the height of the peak from the baseline (y0), x0 is the peak
	center, and w is the standard deviation.  The initial parameter guess are:
		x0, y-weighted mean of the x values,
		sig, sqrt(8*log(2)) * the y-weighted standard deviation of the x values,
		A, maximum y value,
		y0, median value of the y data.
	This function also returns the R**2 "goodness" of the fit (Rsqr).
	
	Args:
		X (list-like): x-axis data
		Y (list-like): y-axis data
	Returns:
		(numpy.ndarray) [x0, w, A, y0, Rsqr]
	'''
	
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	pInit = np.array([0.0, 0.0, max(Y), np.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = np.sqrt(8.0*np.log(2.0)) * np.sqrt( ex1 - pInit[0]**2 )
	
	def lntzResid(params):
		x0 = params[0]
		w = params[1]
		A = params[2]
		y0 = params[3]
		F = A / (1.0 + 4.0*((X-x0)/w)**2 ) + y0
		return F-Y
	
	(fitParams, cov_x, infodict, mesg, intFlag) = leastsq(
		lntzResid, pInit, full_output=True
		)
	if fitParams[1] < 0:
		fitParams[1] = abs(fitParams[1])
	
	return fitParams		
# END lorentzFit

#===============================================================================
def linfit(X, Y, eY=None, full_output=False):
	'''Weighted Linear Best-Fit
	
	Functional form of the linear fit: y(x) = a*x + b
	
	Args:
		X (list-like): x-axis data
		Y (list-like): y-axis data
		eY (list|float|int): Optional, list or value representing the error of
							  each y value
		full_output (bool): If True, extra output will be given
	Returns:
		(tuple) For full_output=False, (a, b).  Otherwise, if error is not
		given, (a, b, Rsqr), where Rsqr is the R**2, else if error is specified,
		(a, b, a_err, b_err, reduc_chisqr), where a_err and b_err are the
		respective fit parameter errors and reduc_chisqr is the reduced chi**2
		value.
	Examples:
		a, b = linfit(X, Y)
		a, b, Rsqr = linfit(X, Y, full_output=True)
		a, b, a_err, b_err, reduc_chisqr = linfit(X, Y, 0.5, True)
		a, b, a_err, b_err, reduc_chisqr = linfit(X, Y, [0.1, 0.2, ...], True)
	'''
	
	if type(X) != list and type(X).__name__ != 'ndarray':
		raise TypeError('X data must be a list or ndarray')
	if type(Y) != list and type(Y).__name__ != 'ndarray':
		raise TypeError('Y data must be a list or ndarray')
	
	n = len(X)
	
	if len(Y) != n:
		raise ValueError('X and Y lists must be the same length')
	
	if eY is None:
		pass
	elif type(eY) == int or type(eY) == float:
		eY = [eY for i in range(0,n)]
	elif type(eY) != list and type(eY).__name__ != 'ndarray':
		raise TypeError('Data error must be an int, float, list, or a ndarray')
	# END if
	
	sx = 0.0
	sx2 = 0.0
	sy = 0.0
	sy2 = 0.0
	sxy = 0.0
	
	if eY is None:
		sw = n
		for i in range(n):
			sx += X[i]
			sx2 += X[i]*X[i]
			sy += Y[i]
			sy2 += Y[i]*Y[i]
			sxy += Y[i]*X[i]
	else:
		for i in range(n):
			sw += 1.0/(eY[i]**2)
			sx += X[i]/(eY[i]**2)
			sx2 += X[i]*X[i]/(eY[i]**2)
			sy += Y[i]/(eY[i]**2)
			sy2 += Y[i]*Y[i]
			sxy += Y[i]*X[i]/(eY[i]**2)
		# END for
	# END if
	
	delta = sw*sx2 - sx*sx
	a = (sw*sxy - sx*sy) / delta
	b = (sx2*sy - sx*sxy) / delta
	
	if not full_output:
		return a, b
	
	if eY is None:
		ssxy = sxy - sx*sy/n
		ssxx = sx2 - sx*sx/n
		ssyy = sy2 - sy*sy/n
		Rsqr = ssxy*ssxy / (ssxx*ssyy)
		return a, b, Rsqr
	else:
		chisqr = 0.0
		for i in range(n):
			chisqr += ( (Y[i] - a*X[i] - b)/eY[i] )**2
		# END for
		reduc_chisqr = chisqr/(n-2)
		
		a_err = sw/delta
		b_err = sx2/delta
		return a, b, a_err, b_err, reduc_chisqr
	# END if
# END linfit
