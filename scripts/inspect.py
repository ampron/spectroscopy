#!/usr/bin/env python2.7
'''Multispec Curves Inspection Script
	
	[description]
	
	Module dependencies:
		matplotlib
		os
		re
		spectroscopy
	List of classes: -none-
	List of functions:
		main
'''

# built-in modules
import os
import re

# third-party modules
import matplotlib.pyplot as plt

# self package
import spectroscopy as spec

#===============================================================================
def main(*args, **kwargs):
	print 'Running "inspect" script'
	
	# Find all multispec (.*vms.asc) files
	specfiles = spec.find_multispec()
	specfiles.sort()
	
	# Curve plotting loop
	type_lkup = {'zv': spec.ZVSTSBundle, 'iv': spec.IVSTSBundle}
	for file_path in specfiles:
		print 'Plotting all from "{}"'.format(file_path)
		path_mat = re.search(r'^(.*?)\.([zi]v)ms(?:\.asc)?$', file_path)
		new_foldername = path_mat.group(1)
		spectype = path_mat.group(2)
		all_curves = type_lkup[spectype](file_path)
		
		# Make new folder for saving individual plots
		if not os.path.exists(new_foldername):
			os.makedirs(new_foldername)
		
		X = all_curves.X
		for i in range(all_curves.N):
			save_path = new_foldername+'/curve_{0:03d}.png'.format(i)
			if os.path.exists(save_path): continue
			
			Y = all_curves.allY[i]
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.plot(X, Y)
			spec.pltformat_basic(ax)
			ax.set_title( 'Curve i={}'.format(i) )
			
			fig.savefig(save_path, dpi=100)
			plt.close(fig)
		# END for
	# END for
	
	print 'inspect is finished.'
# END main

#===============================================================================
if __name__ == '__main__':
	main()
# END of module
