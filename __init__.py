'''AMP Spectroscopy Analysis Package Initialization Module
	Author: Alex Pronschinske
	
	Module dependencies: 
		bundles
		specFit
		specPlot
		specTools
'''

# All functions and classes from the following modules should be imported as
#  native parts of the spectroscopy package
from specFit import *
from specPlot import *
from bundles import *
from specTools import *
from numerical import *

__version__ = '02.00rc'