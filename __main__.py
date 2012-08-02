#!/usr/bin python2.7
'''Spectroscopy Analysis Package Executable Commands Module
	Author: Alex Pronschinske
	
	This module allows the spectroscopy package to be run as script.
	
	List of Classes: -none-
	List of Functions: 
		main
'''

# built-in modules
import sys
import traceback
import re

# internal modules
from spectroscopy import scripts
from spectroscopy.scripts import *

#===============================================================================
def main(*args):
	args = list(args)
	try:
		script_name = args.pop(0)
	except IndexError:
		print 'calling show_menu again'
		args = show_menu()
		main(*args)
		return None
	# END try
	
	if script_name == 'quit': return None
	
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
def show_menu():
	options = scripts.__all__ + ['quit']
	print 'Please enter the name of one of the following scripts or commands:'
	for opt in options:
		print 4*' ' + opt
	# END for
	
	while True:
		user_input = raw_input('--> ')
		user_input = re.split(r'(?<!\\) ', user_input)
		if user_input[0] in options: break
		else: print 'Invalid option'
	# END while
	
	return user_input
# END show_menu

#===============================================================================
if __name__ == '__main__':
	# The following try-except statement is a work-around to allow Windows users
	# to see output and errors when they run the script from the file browser
	# with a double-click
	try:
		main(*sys.argv[1:])
	except Exception as err:
		exc_type, exc_value, exc_tb = sys.exc_info()
		bad_file, bad_line, func_name, text = traceback.extract_tb(exc_tb)[-1]
		print 'Error in {}'.format(bad_file)
		print '{} on {}: {}'.format(type(err).__name__, bad_line, err)
		print ''
	finally:
		raw_input("press any key to exit")
	# END try
# END if