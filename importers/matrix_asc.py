'''MATRIX ASCII File Importing Module
	Author: Alex M. Pronschinske
	
	List of classes: -none-
	List of functions: 
		import_matrix_asc_iv
		import_matrix_asc_zv
'''

# built-in modules
import re

# third-party modules
import numpy as np

#===============================================================================
def import_matrix_asc_iv(file_name):
	'''Function for importing I(V)-STS data
	
	'''
	return _import_matrix_asc(file_name)
# END import_iv

#===============================================================================
def import_matrix_asc_zv(file_name):
	'''Function for importing z(V)-STS data
	'''
	X, allY, units = _import_matrix_asc(file_name)
	
	# conversions from key to nanometers
	to_nm = {'nm': 1.0, 'V': 9.0, 'mV': 9.0E-3, 'nV': 9.0E-9}
	for i in range(len(allY)):
		# scale by 9 nm/V
		allY[i] = allY[i] * to_nm[units[1]]
	# END for
	units[1] = 'nm'
	
	return X, allY, units
# END import_zv

#===============================================================================
def _import_matrix_asc(file_name):
	'''General STS importer function
	'''
	f = open(file_name)
	
	X = []
	for ln in f:
		if re.match(r'[#\r\n]', ln[0]):
			# skip header lines and blank lines
			continue
		elif re.match(r'[a-z]', ln[0], re.IGNORECASE):
			# record units
			colUnits = re.split(r'\s+', ln)
			if colUnits[0] != 'X':
				N = len(colUnits)/2
				allY = [[] for i in range(0,N)]
				units = [colUnits[0], colUnits[1]]
			# END if
		else:
			# record data
			lnData = re.split(r'\s+', ln)
			X.append( float(lnData[0]) )
			for i in range(0,N):
				allY[i].append( float(lnData[2*i+1]) )
		# END if
	# END for loop over file lines
	f.close()
	
	# Reverse the X list if necessary so that it is always in the
	# positive direction
	if X[1]-X[0] < 0: X.reverse()
	
	# convert lists to ndarrays
	X = np.array(X)
	for i in range(len(allY)):
		allY[i] = np.array(allY[i])
		
	return X, allY, units
# END import_file_matrix