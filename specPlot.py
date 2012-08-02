'''Spectroscopy Package Plotting Module
	Author: Alex M. Pronschinske
	
	List of classes: -none-
	List of functions:
		format_default
		plot_density
		tight_scale
'''

# third-party modules
import matplotlib
import matplotlib.pyplot as plt

#===============================================================================
def pltformat_basic(ax):
	'''Format plot: basic style
		
		Formatting adjustments include: axes thickness, axes z-position,
		major & minor tick dimensions, grid appearance, and application of
		matplotlib.tight_layout() function.
		
		Args:
			ax (Axes): matplotlib.Axis object for the plot to be formated
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
	'''DEPRECIATED, use pltformat_basic instead'''
	print (
		'Warning "format_default" will be depreciated, '
		+ 'use "pltformat_basic" instead'
	)
	pltformat_basic(ax)
# END format_default

#===============================================================================
def plot_density(bundle, ax=None, mkrcolor='blue', alpha=None, **extra_kwargs):
	'''Create curve density plot
	
	This function will plot all of the curves in the input bundle on the same
	axis where each is semi-transparent so that where curves overlap the color
	appears darker
	
	Args:
		bundle (SpecBundle): Group of data to be plotted
		ax = None (Axis): Axis on which to add the density plot, if none is
						   given then one will be made and returned
		mkrcolor = 'blue' (str): curve color
		alpha = None (float): The Opacity value of each curve, if none is given
							   if will be set to 0.1 for bundle of 100 or fewer
							   curves and set up 8/bundle.N if their are more
							   than 100 curves
		**extra_kwargs: Any extra keywords will be passed onto the ax.plot()
						function call
	Returns:
		(Axis) The newly created plot
	Example:
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plot_density(bun, ax)
		ax2 = plot_density(bun)
	'''
	
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
	'''Tighten the y-axis scaling
	
	This function will adjust the y-axis limits of a plot so that the specified
	fit_amount fraction will be within the limits then grow the size of the
	vertical frame by the growth_factor fractional amount.  For example the
	default values of fit_amount=0.99 and growth_factor=0.20 will first set the
	y-limits to accommodate 99% of the data, then make the frame +20% bigger
	(i.e. 120% of it's adjusted size).
	
	Args:
		ax (Axis): Axis which will be re-scaled
		fit_amount = 0.99 (float):
		growth_factor = 0.20 (float):
	Returns:
		(tuple) (low_bnd, upp_bnd) New upper and lower bounds of the y-axis.
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