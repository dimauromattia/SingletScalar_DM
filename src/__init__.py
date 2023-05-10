'''
Singlet Scalar DM
=================

This repository should be installed through `pip`, running the following command
in the root directory.

.. code-block::shell
    python -m pip install .

Importing the package automatically loads a matplotlibrc file defining the
plotting style.

Functions
---------


Variables
---------
fN
mN
mh
GeVm2tocm2
Gamma_H_SM
v
mass_vector_micro
mass_vector_drake
massz_vec
MassDD_vec
logenergyx_bins
lambda_vec
lambdahs_vec
'''

from . import _globals
from ._globals import *
from . import _functions
from ._functions import *

__all__ = []
__all__ += _globals.__all__
__all__ += _functions.__all__
