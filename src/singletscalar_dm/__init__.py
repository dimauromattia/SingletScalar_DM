'''
Singlet Scalar DM
=================

This repository should be installed through `pip`, running the following command
in the root directory.

.. code-block:: bash

    python -m pip install .

Functions
---------
interpolate_Omega
interpolate_Omega_MicrOMEGAs
interpolate_lambda
interpolate_lambda_MicrOMEGAs
sigmav_channels
DMspectra_inttable
provide_ULEXP
SI_noomega
SI_withomega
GetUL_DD_nomega
GetUL_DD_withomega
minimize_br_inv
Gamma_inv
Br_inv
Br_inv_UL

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

Notes
-----
The package routines rely on additional data files containing the main data
the computations are based on.
These data files can be obtained by using `importlib.resources`:

>>> from importlib.resources import files
... data_text = files('singletscalar_dm.data').joinpath('SHP_sigmav_bb.dat')
'''

from . import _globals
from ._globals import *
from . import _functions
from ._functions import *

__all__ = []
__all__ += _globals.__all__
__all__ += _functions.__all__
