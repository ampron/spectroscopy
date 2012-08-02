#!/usr/bin/env python2.7
'''STS Quick Batch Plotting Script
	
	Basic script for batch processing of data into "quick plots"
	
	Module dependencies:
		matplotlib
		multiprocessing
		os
		re
		spectroscopy
		sys
	List of classes: -none-
	List of functions:
		main
		zv_subroutine
		iv_subroutine
'''

# built-in modules
import sys
import traceback
import os
import re
import multiprocessing as mproc

# third-party modules
import matplotlib.pyplot as plt

# self package
import spectroscopy as spec

#===============================================================================
def main(path='./', flags=set()):
	print sys.version
	print 'Running "quick_plots" script in {}'.format(path)
	
	# Handle recursive search option
	recursive = False
	if 'r' in flags:
		recursive = True
		print '*Including all sub-folders'
	# END if
	
	# Handle the no-rewrite option
	only_new = False
	if 'n' in flags:
		only_new = True
		print '*Making new plots only'
	# END if
	
	# Handle the multiprocessing options
	top_gear = False
	first_gear = False
	if 'm' in flags:
		first_gear = True
		print '*Optioned to run everything under a single process'
	elif 'M' in flags:
		top_gear = True
		print '*Optioned to run full number of processes'
	# END if
	
	# Find all multispec (.*vms.asc) files
	specfiles = spec.find_multispec(path, recursive)
	specfiles.sort()
	
	# For no-rewrite option create a list of exclusions
	if only_new:
		skips = set()
		files_stack = ['./'+f for f in os.listdir(path)]
		while len(files_stack) > 0:
			target = files_stack.pop()
			name_match = re.search(
				r'^.*?([^/]+?)-(?:zv|dzdv|iv|didv|ndidv)\.png$', target
			)
			if recursive and os.path.isdir(target):
				for f in os.listdir(target): files_stack.append(target+'/'+f)
			elif name_match and name_match.group(1) not in skips:
				skips.add( name_match.group(1) )
			# END if
		# END while
	# END if
	
	# Setup the worker processes pool
	N_procs = mproc.cpu_count()
	if first_gear:
		N_procs = 1
	elif not top_gear:
		# Default to using one less process than the total number of processors,
		# for safety
		N_procs -= 1
		if N_procs == 0:
			N_procs = 1
			print (
				'*Default safe running conditions will only use'
				+ 'a single process on this system'
			)
		# END if
	# END if
	pool = mproc.Pool(processes=N_procs, maxtasksperchild=12)
	
	# Plotting loop
	for file_path in specfiles:
		path_mat = re.search(r'^(.+?)([^/]+)$', file_path)
		path = path_mat.group(1)
		file_name = path_mat.group(2)
		
		if (only_new
		and re.sub(r'\.[zi]vms(?:\.asc)?', '', file_name) in skips):
			continue
		
		type_mat = re.search(r'\.([zi]v)ms(?:\.asc)?$', file_name)
		spectype = type_mat.group(1)
		if spectype == 'zv':
			pool.apply_async(zv_subroutine, args=(path, file_name))
		elif spectype == 'iv':
			if os.path.exists(re.sub(r'\.ivms', '.gvms', file_path)):
				pool.apply_async(giv_subroutine, args=(path, file_name))
			else:
				pool.apply_async(iv_subroutine, args=(path, file_name))
		else:
			print "don't know what to do with {}".format(file_name)
		# END if
	# END for
	
	pool.close()
	pool.join()
	
	print 'quick_plots is finished.'
# END main

#===============================================================================
def zv_subroutine(path, file_name):
	try:
		_zv_subroutine(path, file_name)
	except Exception as err:
		exc_type, exc_value, exc_tb = sys.exc_info()
		bad_file, bad_line, func_name, text = traceback.extract_tb(exc_tb)[-1]
		print 'Error processing {}'.format(path+file_name)
		print '  from {}'.format(bad_file)
		print '  {} on {}: {}'.format(type(err).__name__, bad_line, err)
		print ''
		raise Exception
	# END try
# END zv_subroutine

#===============================================================================
def _zv_subroutine(path, file_name):
	log = 'z(V) spectra "{}":\n'.format(file_name)
	log += '  folder: {}\n'.format(path)
	
	allZs = spec.ZVSTSBundle(path+file_name)
	
	rr = allZs.has_ramp_reversal
	if rr:
		log += '  ramp reversal detected\n'
		allZs_f, allZs_r = allZs.ramp_reversal_fix(opt='split')
		log += '  keeping track of ramp directions\n'
		
		N_spec = allZs_f.N
		mZ_f = allZs_f.coavg()
		mZ_r = allZs_r.coavg()
		V = allZs_f.X
	else:
		N_spec = allZs.N
		mZ = allZs.coavg()
		V = allZs.X
	# END if
	
	log += '  {} spectra in sample\n'.format(N_spec)
	
	# make density plot for all the z(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(allZs_f, ax)
		spec.plot_density(allZs_r, ax, mkrcolor='red')
		ax.plot(V, mZ_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, mZ_r, color=(1.0,0.5,0.5), linewidth=2)
	else:
		spec.plot_density(allZs, ax)
		ax.plot(V, mZ, color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	spec.pltformat_basic(ax)
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-zv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	# create a dz/dV(V) spectra bundle
	x_window = 0.1 #V
	poly_order = 1
	deriv_order = 1
	if rr:
		alldZs_f = allZs_f.deriv_sg(x_window, poly_order, deriv_order)
		alldZs_r = allZs_r.deriv_sg(x_window, poly_order, deriv_order)
		mdZ_f = alldZs_f.coavg()
		mdZ_r = alldZs_r.coavg()
	else:
		alldZs = allZs.deriv_sg(x_window, poly_order, deriv_order)
		mdZ = alldZs.coavg()
	# END if
	
	# make density plot for all the dz/dV(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(alldZs_f, ax)
		spec.plot_density(alldZs_r, ax, mkrcolor='red')
		ax.plot(V, mdZ_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, mdZ_r, color=(1.0,0.5,0.5), linewidth=2)
	else:
		spec.plot_density(alldZs, ax)
		ax.plot(V, mdZ, color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	spec.pltformat_basic(ax)
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-dzdv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	print log
# END _zv_subroutine

#===============================================================================
def iv_subroutine(path, file_name):
	try:
		_iv_subroutine(path, file_name)
	except Exception as err:
		exc_type, exc_value, exc_tb = sys.exc_info()
		bad_file, bad_line, func_name, text = traceback.extract_tb(exc_tb)[-1]
		print 'Error processing {}'.format(path+file_name)
		print '  from {}'.format(bad_file)
		print '  {} on {}: {}'.format(type(err).__name__, bad_line, err)
		print ''
		raise Exception
	# END try
# END iv_subroutine

#===============================================================================
def _iv_subroutine(path, file_name):
	log = 'I(V) spectra "{}":\n'.format(file_name)
	log += '  folder: {}\n'.format(path)
	
	allIs = spec.IVSTSBundle(path+file_name)
	
	rr = allIs.has_ramp_reversal
	if rr:
		log += '  ramp reversal detected\n'
		allIs_f, allIs_r = allIs.ramp_reversal_fix(opt='split')
		log += '  keeping track of ramp directions\n'
		
		N_spec = allIs_f.N
		mI_f = allIs_f.coavg()
		mI_r = allIs_r.coavg()
		V = allIs_f.X
	else:
		N_spec = allIs.N
		mI = allIs.coavg()
		V = allIs.X
	# END if
	log += '  {} spectra in sample\n'.format(N_spec)
	
	# make density plot for all the I(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(allIs_f, ax)
		spec.plot_density(allIs_r, ax, mkrcolor='red')
		ax.plot(V, mI_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, mI_r, color=(1.0,0.5,0.5), linewidth=2)
	else:
		spec.plot_density(allIs, ax)
		ax.plot(V, mI, color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	spec.pltformat_basic(ax)
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-iv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	# create a dz/dV(V) spectra bundle
	x_window = 0.1 #V
	poly_order = 1
	if rr:
		allndIs_f = allIs_f.norm_deriv(x_window, poly_order)
		allndIs_r = allIs_r.norm_deriv(x_window, poly_order)
		dv = V[1] - V[0]
		ndmI_f = spec.norm_deriv(V, allIs_f.coavg(), x_window, poly_order)
		ndmI_r = spec.norm_deriv(V, allIs_r.coavg(), x_window, poly_order)
	else:
		allndIs = allIs.norm_deriv(x_window, poly_order)
	# END if
	
	# make density plot for all the norm. dI/dV(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(allndIs_f, ax)
		spec.plot_density(allndIs_r, ax, mkrcolor='red')
		ax.plot(V, ndmI_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, ndmI_r, color=(1.0,0.5,0.5), linewidth=2)
		ax.plot(V, allndIs_f.coavg(), color=(0.5,0.5,0.6), linewidth=1)
		ax.plot(V, allndIs_r.coavg(), color=(0.6,0.5,0.5), linewidth=1)
	else:
		spec.plot_density(allndIs, ax)
		ax.plot(V, allndIs.coavg(), color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	spec.pltformat_basic(ax)
	spec.tight_scale(ax)
	
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-ndidv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	print log
# END _iv_subroutine

#===============================================================================
def giv_subroutine(path, file_name):
	try:
		_giv_subroutine(path, file_name)
	except Exception as err:
		exc_type, exc_value, exc_tb = sys.exc_info()
		bad_file, bad_line, func_name, text = traceback.extract_tb(exc_tb)[-1]
		print 'Error processing {}'.format(path+file_name)
		print '  from {}'.format(bad_file)
		print '  {} on {}: {}'.format(type(err).__name__, bad_line, err)
		print ''
		raise Exception
	# END try
# END giv_subroutine

#===============================================================================
def _giv_subroutine(path, file_name):
	log = 'I(V) spectra "{}":\n'.format(file_name)
	gspec_file_name = re.sub(r'\.ivms', '.gvms', file_name)
	log += '  dI/dV(V) spectra: "{}"'.format(gspec_file_name)
	log += '  folder: {}\n'.format(path)
	
	allIs = spec.IVSTSBundle(path+file_name)
	alldIs = spec.IVSTSBundle(path+gspec_file_name)
	
	rr = allIs.has_ramp_reversal
	if rr:
		log += '  ramp reversal detected\n'
		allIs_f, allIs_r = allIs.ramp_reversal_fix(opt='split')
		alldIs_f, alldIs_r = alldIs.ramp_reversal_fix(opt='split')
		log += '  keeping track of ramp directions\n'
		
		N_spec = allIs_f.N
		mI_f = allIs_f.coavg()
		mI_r = allIs_r.coavg()
		mdI_f = alldIs_f.coavg()
		mdI_r = alldIs_r.coavg()
		V = allIs_f.X
	else:
		N_spec = allIs.N
		mI = allIs.coavg()
		mdI = alldIs.coavg()
		V = allIs.X
	# END if
	log += '  {} spectra in sample\n'.format(N_spec)
	
	# make density plot for all the I(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(allIs_f, ax)
		spec.plot_density(allIs_r, ax, mkrcolor='red')
		ax.plot(V, mI_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, mI_r, color=(1.0,0.5,0.5), linewidth=2)
	else:
		spec.plot_density(allIs, ax)
		ax.plot(V, mI, color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	ax.set_xlabel('Gap Bias ({})'.format(allIs.units[0]))
	ax.set_ylabel('Tunneling Current ({})'.format(allIs.units[1]))
	spec.pltformat_basic(ax)
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-iv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	# make density plot for all the dI/dV(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(alldIs_f, ax)
		spec.plot_density(alldIs_r, ax, mkrcolor='red')
		ax.plot(V, mdI_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, mdI_r, color=(1.0,0.5,0.5), linewidth=2)
	else:
		spec.plot_density(alldIs, ax)
		ax.plot(V, mdI, color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	ax.set_xlabel('Gap Bias ({})'.format(alldIs.units[0]))
	ax.set_ylabel('Lock-in Signal ({})'.format(alldIs.units[1]))
	spec.pltformat_basic(ax)
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-didv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	# create a dz/dV(V) spectra bundle
	x_window = 0.1 #V
	poly_order = 1
	if rr:
		sIs_f = allIs_f.deriv_sg(x_window, poly_order, 0)
		sIs_f.zero_y0()
		sIs_r = allIs_r.deriv_sg(x_window, poly_order, 0)
		sIs_r.zero_y0()
		sdIs_f = alldIs_f.deriv_sg(x_window, poly_order, 0)
		sdIs_r = alldIs_r.deriv_sg(x_window, poly_order, 0)
		allndIs_f = alldIs_f.copy()
		allndIs_r = alldIs_r.copy()
		for i in range(allndIs_f.N):
			allndIs_f[i] = allIs_f.X * sdIs_f[i] / sIs_f[i]
			allndIs_r[i] = allIs_f.X * sdIs_r[i] / sIs_r[i]
		# END for
		
		ndmI_f = spec.norm_deriv(V, allIs_f.coavg(), x_window, poly_order)
		ndmI_r = spec.norm_deriv(V, allIs_r.coavg(), x_window, poly_order)
	else:
		sIs = allIs.deriv_sg(x_window, poly_order, 0)
		sIs.zero_y0()
		sdIs = alldIs.deriv_sg(x_window, poly_order, 0)
		allndIs = alldIs.copy()
		for i in range(allndIs.N):
			allndIs[i] = allIs.X * sdIs[i] / sIs[i]
		# END for
	# END if
	
	# make density plot for all the norm. dI/dV(V) spectra
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if rr:
		spec.plot_density(allndIs_f, ax)
		spec.plot_density(allndIs_r, ax, mkrcolor='red')
		ax.plot(V, ndmI_f, color=(0.5,0.6,1.0), linewidth=2)
		ax.plot(V, ndmI_r, color=(1.0,0.5,0.5), linewidth=2)
		ax.plot(V, allndIs_f.coavg(), color=(0.5,0.5,0.6), linewidth=1)
		ax.plot(V, allndIs_r.coavg(), color=(0.6,0.5,0.5), linewidth=1)
	else:
		spec.plot_density(allndIs, ax)
		ax.plot(V, allndIs.coavg(), color=(0.5,0.6,1.0), linewidth=2)
	# END if
	ax.set_title( '{}, N={}'.format(file_name, N_spec) )
	ax.set_xlabel('Gap Bias ({})'.format(alldIs.units[0]))
	ax.set_ylabel('(V/I)*dI/dV')
	spec.pltformat_basic(ax)
	spec.tight_scale(ax)
	
	save_name = re.sub(r'\.[zi]vms(?:\.asc)?$', '-ndidv.png', file_name)
	ax.figure.savefig(path+save_name, dpi=150)
	plt.close(fig)
	log += '  saved "{}"\n'.format(save_name)
	
	print log
# END _giv_subroutine

#===============================================================================
if __name__ == '__main__':
	main(*sys.argv[1:])
# END of module
