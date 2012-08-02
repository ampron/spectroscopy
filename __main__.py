#!/usr/bin/env python2.7
'''AMP Spectroscopy Analysis Package Executable Commands Module
	Author: Alex Pronschinske
	
	List of Classes: -none-
	List of Functions: 
		main
	Module dependencies: 
		spectroscopy
		sys
'''

import sys
from spectroscopy.scripts import *

#===============================================================================
def main(*args):
	args = list(args)
	script_name = args.pop(0)
	
	flags = set()
	if len(args) > 0 and args[0][0] == '-':
		flags.update( *list(args.pop(0)[1:]) )
	
	if script_name in globals():
		globals()[script_name].main(*args, flags=flags)
	else:
		print 'script "{}" does not exist'.format(script_name)
	# END if
# END main

#===============================================================================
if __name__ == '__main__':
	main(*sys.argv[1:])