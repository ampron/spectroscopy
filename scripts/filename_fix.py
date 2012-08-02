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
import os

# self package
import spectroscopy as spec

#===============================================================================
def main_alt(path='./', flags=set()):
	print 'Running "filename_fix" script'
	
	# Simulate flag
	simulate = False
	if 's' in flags:
		simulate = True
		print '*Simulating file name changes'
	# END if
	
	# Find all multispec (.*vms.asc or .*vms) files
	specfiles = spec.find_multispec(path)
	specfiles.sort()
	
	# file renaming loop
	for file_path in allfiles:
		fname = os.path.basename(file_path)
		print '"{}"'.format(fname)
		
		tags = re.search(r'''
			^( \d{6}                                         ) #0 date
			(  (?:-\d{4})?                                   ) #1 time
			(  (?:-[^-]+)*?                                  ) #2 pre-index
			(  -\d+ (?:_\d+ (?:_[fubd]+)?)? (?:-[A-Z]{1,2})? ) #3 index
			(  (?: -[^-]+ )*                                 ) #4 post-index
			(  \.[zi]vms (?:\.asc)?                         )?$ #5 file-extention
			''', fname, re.VERBOSE
		).groups()
		
		new_name = tags[0] + tags[1] + tags[2] + tags[4]
		new_name += re.sub(r'(?<!^)-', '_', tags[3])
		new_name += re.sub(r'\.asc', '', tags[5])
		
		print '  -->"{}"'.format(new_name)
		if not simulate and new_name != fname: os.rename(fname, new_name)
	# END for
	
	print 'filename_fix is finished.'
# END main

#===============================================================================
if __name__ == '__main__':
	main()
# END of module
