'''AMP Data Fitting Module
	Author: Alex M. Pronschinske
	
	List of classes: 
		PeakFitter
	List of functions:
		gaussFit
		gaussLinFit
		lorentzFit
		linfit
	Module dependencies: 
		scipy
		scipy.optimize
'''

import scipy as sp
from scipy.optimize import leastsq

class PeakFitter():
	'''Peak fitting classes
	
	Class methods:
	'''
	
	def __init__():
		pass
	# END __init__
	
	def guess(self, X, Y):
		'''Peak Parameter Guessing Method
		
		Args:
			X (list-like): X data points
			Y (list-like): Y data points
		Returns:
			(tuple)  (x0_guess, std_guess)
		'''
		
		Y = (Y-min(Y)) / sum(Y)
		x0_guess = 0.0
		std_guess = 0.0
		sum_sqrs = 0.0
		for i in range(0, len(X)):
			x0_guess += X[i]*Y[i]
			sum_sqrs += X[i]**2 * Y[i]
		std_guess = sp.sqrt( sum_sqrs - x0_guess**2 )
		
		return (x0_guess, std_guess)
	# END guess
# END PeakFitter

def gaussFit(X, Y):
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	pInit = sp.array([0.0, 0.0, max(Y), sp.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = sp.sqrt( ex1 - pInit[0]**2 )
	
	def GaussResid(params):
		x0 = params[0]
		sig = params[1]
		A = params[2]
		y0 = params[3]
		F = A * sp.exp(-0.5*((X-x0)/sig)**2 ) + y0
		return F-Y
	
	(fitParams, cov_x, infodict, mesg, intFlag) = leastsq(
		GaussResid, pInit, full_output=True
		)
	if fitParams[1] < 0:
		fitParams[1] = abs(fitParams[1])
	
	Yf = infodict['fvec']
	Rsq = 1.0 - ( sp.std(Y)**2/sum((Y-Yf)**2) )
	
	return sp.append(fitParams, Rsq)
# END gaussFit

def gaussLinFit(X, Y):
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	#Parameter Array: [x0, sig, A, y1, y0]
	pInit = sp.array([0.0, 0.0, max(Y), 0.0, sp.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = sp.sqrt( ex1 - pInit[0]**2 )
	
	def GaussResid(params):
		x0 = params[0]
		sig = params[1]
		A = params[2]
		y1 = params[3]
		y0 = params[4]
		F = A * sp.exp(-0.5*((X-x0)/sig)**2 ) + y1*X + y0
		return F-Y
	
	(fitParams, cov_x, infodict, mesg, intFlag) = leastsq(
		GaussResid, pInit, full_output=True
		)
	if fitParams[1] < 0:
		fitParams[1] = abs(fitParams[1])
	
	Yf = infodict['fvec']
	Rsq = 1.0 - ( sp.std(Y)**2/sum((Y-Yf)**2) )
	
	return sp.append(fitParams, Rsq)
# END gaussFit

def lorentzFit(X, Y):
	# Find initial guess for x0 and w
	Y2 = Y - min(Y)
	Y2 = Y2 / sum(Y2)
	pInit = sp.array([0.0, 0.0, max(Y), sp.median(Y)])
	ex1 = 0.0
	for i in range(0, len(X)):
		pInit[0] = pInit[0] + X[i]*Y2[i]
		ex1 += X[i]**2 * Y2[i]
	pInit[1] = sp.sqrt(8.0*sp.log(2.0)) * sp.sqrt( ex1 - pInit[0]**2 )
	
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

def linfit(X, Y, eY=None, full_output=False):
	'''Weighted Linear Best-Fit
	
	Args:
		X (list):
		Y (list):
		eY (list|float|int):
		full_output (bool):
	Returns:
		
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
		return (a, b)
	
	if eY is None:
		ssxy = sxy - sx*sy/n
		ssxx = sx2 - sx*sx/n
		ssyy = sy2 - sy*sy/n
		Rsqr = ssxy*ssxy / (ssxx*ssyy)
		return (a, b, Rsqr)
	else:
		chisqr = 0.0
		for i in range(n):
			chisqr += ( (Y[i] - a*X[i] - b)/eY[i] )**2
		# END for
		reduc_chisqr = chisqr/(n-2)
		
		a_err = sw/delta
		b_err = sx2/delta
		return (a, b, a_err, b_err, reduc_chisqr)
	# END if
# END linfit
