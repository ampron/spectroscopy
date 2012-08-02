'''Importers Sub-Package Initialization Module
	Author: Alex Pronschinske
	
	List of Classes: -none-
	List of Functions: -none-
'''

# Built-in modules
import sys
import os

# This will setup the __all__ list dynamically, so that it includes all modules
# in this sub-package.
# This allows for "from spectroscopy.importers import *" to import all script
# modules that are in the scripts folder
__all__ = []
for file_name in os.listdir(os.path.dirname(__file__)):
	if file_name[-3:] == '.py' and file_name != '__init__.py':
		__all__.append(file_name[:-3])
# END for
del file_name

# This allows for "from spectroscopy import importers" to import an importers
# module that will contain all of the modules in the importers folder
sys.path.append(os.path.dirname(__file__))
for mod_name in __all__:
	globals()[mod_name] = __import__(mod_name)
del mod_name