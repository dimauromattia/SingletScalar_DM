'''
Singlet Scalar DM
=================

This modules provides the code for calculating cross sections, relic density
and source spectra for the singlet scalar model of dark matter.

Notes
-----
The package routines rely on additional data files containing the main data
the computations are based on.
These data files can be obtained by using the function `import_data_file`:

>>> from singletscalar_dm import import_data_file
... data_path = import_data_file('SHP_sigmav_bb.dat')

References
----------


Examples
--------
The ``examples`` folder contains python examples ready to be run to plot relevant
phenomenological quantities, leveraging the functions defined in this module.
The online documentation provides an handy way of browsing these examples,
enabling the possibility of downloading them as `.py` scripts or `.ipynb` notebooks.
'''

from . import _globals
from ._globals import *
from . import _functions
from ._functions import *

__all__ = []
__all__ += _globals.__all__
__all__ += _functions.__all__
