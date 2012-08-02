#!/usr/bin/env python2.7
'''Multispec File Merging Script
	
	[description]
	
	Module dependencies:
		matplotlib
		multiprocessing
		os
		re
		spectroscopy
	List of classes: -none-
	List of functions:
		main
'''

# built-in modules
import re

# self package
import spectroscopy as spec

#===============================================================================
def main(*args, **kwargs):
	print 'Running "ms_merge" script'
	
	# Find all multispec (.*vms.asc) files
	specfiles = spec.find_multispec()
	specfiles.sort()
	
	# file reading and validating loop
	exp_spectype = None
	new_bundle = None
	type_lkup = {'zv': spec.ZVSTSBundle, 'iv': spec.IVSTSBundle}
	for file_path in specfiles:
		try:
			spectype = re.search(r'\.([zi]v)ms(?:\.asc)?$', file_path).group(1)
			if exp_spectype is None:
				exp_spectype = spectype
			elif exp_spectype != spectype:
				print 'Cannot mix spectroscopy types'
				print 'ms_merge failed.'
				return None
			# END if
		except AttributeError:
			print '"{}" is not a multispec file'.format(file_path)
			print 'ms_merge failed.'
			return None
		# END try
		
		ms_data = type_lkup[spectype](file_path)
		
		if new_bundle is None:
			new_bundle = ms_data
		elif list(new_bundle.X) != list(ms_data.X):
			print 'Cannot concatenate multispecs with different domains'
			print '{}:'.format(file_path)
			print '  x_0 = {}'.format(ms_data.X[0])
			print '  x_{} = {}'.format(ms_data.curve_len, ms_data.X[-1])
			print 'Expected:'.format(file_path)
			print '  x_0 = {}'.format(new_bundle.X[0])
			print '  x_{} = {}'.format(new_bundle.curve_len, new_bundle.X[-1])
			print 'ms_merge failed.'
			return None
		# END if
		
		print u'\u2714 {}'.format(file_path)
		new_bundle.concat(ms_data)
	# END for
	
	# Create file name
	# ask for file name from user
	new_name = raw_input('Please enter name for new multispec file:\n-->')
	# remove any path syntax, force save in '.'
	new_name = re.sub(r'^(.*?)([^/]+)$', r'\2', new_name)
	# ensure the correct file extention, ".[zi]vms"
	if not re.search(r'\.[zi]vms$', new_name):
		new_name += '.'+exp_spectype+'ms'
	# END if
	new_bundle.save_multispec(new_name)
	print 'saved as "{}"'.format(new_name)
	
	print 'ms_merge is finished.'
# END main

#===============================================================================
if __name__ == '__main__':
	main()
# END of module
