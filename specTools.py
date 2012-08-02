'''Spectroscopy Package Tools Module
	Author: Alex Pronschinske
	
	List of Classes: -none-
	List of Functions:
		find_multispec
		spipNameFix
		split_spectrum
		re
'''

# built-in modules
import os
import re

# third-party modules
import numpy as np

#===============================================================================
def find_multispec(directory='./', r=False, spec_filter=None):
	'''Find SPIP-formated MATRIX ASCII multiple spectroscopy ("multispec") files
	
	Args:
		directory (str): string of path location
		r (bool): "r"ecursive search, if True include all sub-folders
		spec_filter (str): 'zv', 'iv', or any vaild regular expression pattern.
							Function will return only files with names that
							match the given pattern.
	Returns:
		(list) ['./file', './file', './folder/file']
	'''
	if directory[-1] != '/': directory += '/'
	
	dir_contents = os.listdir(directory)
	
	if spec_filter is None:
		pattern = r'^.+?\.[iz]vms(?:\.asc)?$'
	elif spec_filter.lower() == 'zv':
		pattern = r'^.+?\.zvms(?:\.asc)?$'
	elif spec_filter.lower() == 'iv':
		pattern = r'^.+?\.ivms(?:\.asc)?$'
	else:
		pattern = spec_filter
	# END
	
	dir_contents.insert(0, None)
	while dir_contents[-1] is not None:
		target = dir_contents.pop()
		# END if
		if r and os.path.isdir(directory+target):
			dir_contents = (
				find_multispec(directory+target+'/', True, spec_filter)
				+ dir_contents
			)
		elif re.search(pattern, target):
			dir_contents.insert(0, directory+target)
		# END if
	# END while
	dir_contents.pop()
	
	return dir_contents
# END find_multispec

#===============================================================================
def spipNameFix(dirName = './'):
	if dirName[-1] != '/':
		dirName += '/'
	
	allFiles = os.listdir(dirName)
	
	months = {
		'jan': '01', 'feb': '02', 'mar': '03', 'apr': '04', 'may': '05',
		'jun': '06', 'jul': '07', 'aug': '08', 'sep': '09', 'oct': '10',
		'nov': '11', 'dec': '12'
	}
	
	scDirKey = {
		'TraceUp': 'fu', 'RetraceUp': 'bu',
		'TraceDown': 'fd', 'RetraceDown': 'bd'
	}
	
	# Example file name that will match...
	# ""
	fNamePtn = re.compile(
		r'''
		z \s
		(\w+) \s             #(1) scan direction
		\w\w\w \s            #day of the week
		(\w\w\w) \s          #(2) month
		(\d{1,2}) \s         #(3) date
		([.\d]+?) \s         #(4) time
		(\d{4}) \s           #(5) year
		\[ (\d+-\d+) \] \s   #(6) index
		.+? \.jpg            #end
		''', re.VERBOSE | re.IGNORECASE
		)
	
	for fName in allFiles:
		fNameMat = fNamePtn.match(fName)
		if fNameMat:
			yy = fNameMat.group(5)[2:]
			mm = months[fNameMat.group(2).lower()]
			dd = '{0:2d}'.format(int(fNameMat.group(3)))
			time = re.sub('\.', '', fNameMat.group(4))[:-2]
			index = re.sub('-', '_', fNameMat.group(6))
			scDir = scDirKey[fNameMat.group(1)]
			newFName = yy + mm + dd + '-' + time + '-' + index + '_' \
				+ scDir + '.jpg'
			os.rename(dirName+fName, dirName+newFName)
# END spipNameFix

#===============================================================================
def split_spectrum(X, Y):
	'''Separate a horizontally concatenated (i.e. ramp-reversed) pair of curves
		
		Args:
			X (list-like): x-axis data
			Y (list-like): y-axis data
		Returns:
			(tuple) (X, lY, rY) New x-axis, y-data from the left half (lY),
			and y-data from the right half (rY).
	'''
	
	if len(X)%2 == 1:
		raise ValueError(
			'Spectrum must length must be even, given len(X)='+str(len(X))
		)
	# END if
	
	h = len(X)/2
	# "l"eft Y array
	lY = Y[:h]
	# "r"ight Y array
	rY = Y[h:]
	
	X = np.linspace(X[0], X[-1], len(X)/2)
	
	return X, lY, rY
# split_spectrum