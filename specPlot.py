'''Custom Plotting Module
	Author: Alex M. Pronschinske
	Version: 1
	
	List of classes: -none-
	List of functions:
		format_default
		plot_density
		tight_scale
	Module dependencies: 
		matplotlib
		matplotlib.pyplot
'''

import matplotlib
import matplotlib.pyplot as plt

#===============================================================================
def pltformat_basic(ax):
	'''Format plot: basic style
		
		Formatting adjustments include: axes thickness, axes z-position,
		major & minor tick dimensions, grid appearance, and application of
		matplotlib.tight_layout() function.
		
		Args:
			ax (Axes): Representative object for the plot to be formated
		Returns:
			N/A
	'''
	ax.set_axisbelow(True)
	
	for child in ax.get_children():
		if isinstance(child, matplotlib.spines.Spine):
			child.set_linewidth(2)
	# END for
	ax.tick_params(which='major', labelsize=18, length=7, width=3)
	ax.tick_params(which='minor', length=5, width=3)
	#ax.set_xlabel(size=22)
	#ax.set_ylabel(size=22)
	ax.grid(color=(0.8,0.8,0.8), linewidth=1.5)
	ax.minorticks_on()
	
	ax.figure.tight_layout()
# END pltformat_basic

#===============================================================================
def format_default(ax):
	print (
		'Warning "format_default" will be depreciated, '
		+ 'use "pltformat_basic" instead'
	)
	pltformat_basic(ax)
# END format_default

#===============================================================================
def plot_density(bundle, ax=None, mkrcolor='blue', alpha=None, **extra_kwargs):
	if ax is None:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	# END if
	
	if mkrcolor == 'blue':
		mkrcolor = (0.5,0.6,1.0)
	elif mkrcolor == 'red':
		mkrcolor = (1.0,0.5,0.5)
	# END if
	
	if alpha is None:
		alpha = 0.1
		if bundle.N >= 100: alpha = 8.0/bundle.N
	# END if
	
	X = bundle.X
	for Y in bundle:
		ax.plot(
			X, Y, '-', color=mkrcolor, alpha=alpha, linewidth=2, markersize=5,
			markeredgewidth=0, **extra_kwargs
		)
	
	return ax
# END plot_density

#===============================================================================
def tight_scale(ax, fit_amount=0.99, growth_factor=0.20):
	'''TODO: need documentation
	'''
	
	# Put all of the y-points of all of the curves into one array and sort them
	all_curves = ax.get_lines()
	ypoints = []
	for curve in all_curves: ypoints.extend( curve.get_ydata() )
	ypoints.sort()
	
	# Scale the y-limits
	oneside_lim = int((1-fit_amount)*len(ypoints))/2
	y_window = abs(ypoints[oneside_lim] - ypoints[-oneside_lim])
	low_bnd = ypoints[oneside_lim] - growth_factor*y_window/2
	upp_bnd = ypoints[-oneside_lim] + growth_factor*y_window/2
	
	# Set the new y-limits
	ax.set_ylim(low_bnd, upp_bnd)
	
	# Return the new y-limits in case the user wants them
	return low_bnd, upp_bnd
# END tight_scale