.. SingletScalar-DM documentation master file, created by
   sphinx-quickstart on Thu May 11 00:07:44 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SingletScalar-DM's documentation!
============================================

Installation
============

The package supports Python 3.7 and beyond.
Moreover, it depends on the following packages:

- |numpy|_
- |scipy|_

To install the package along its dependencies, download a release (for stable
version) or clone this repository (for development version) and run the following
command in the package directory:

.. code-block:: bash

    python -m pip install .

.. |numpy| replace:: ``numpy``
.. _numpy: https://numpy.org/
.. |scipy| replace:: ``scipy``
.. _scipy: https://scipy.org/

Running the examples
--------------------

We provide an examples gallery, where it is possible to see many applications
of the routines contained in the package to plot relevant phenomenological
quantities.
It also enables the possibility of downloading every example either as `.py`
scripts or `.ipynb` notebooks.

The examples rely on the |matplotlib|_ package in order to show the plots.
You can install this dependence using the following command in the package root
directory:

.. code-block:: bash

    python -m pip install '.[examples]'

.. |matplotlib| replace:: ``matplotlib``
.. _matplotlib: https://matplotlib.org/

Building the documentation
--------------------------

If you want to rebuild the documentation for yourself, you need to use |sphinx|_
with the configuration provided in the `docs/` folder.

You can install the required dependencies using the following command in the package
root directory:

.. code-block:: bash

    python -m pip install '.[docs]'

Build as usual following the |sphinx|_ directives.

.. |sphinx| replace:: ``sphinx``
.. _sphinx: https://www.sphinx-doc.org/en/master/

Contents
========

.. toctree::
   Back to Github <https://github.com/dimauromattia/SingletScalar_DM>

.. toctree::
   :maxdepth: 1
   :caption: Topics:

   modules
   singletscalar_dm
   examples_gallery/index

Citation
========

If you are using this tool or the results obtained with it in a published work,
please cite `arXiv:####.##### <https://arxiv.org/abs/####.#####>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
==========

The package relies on pre-computed data obtained according to several codes:

- MadDM [MadDM1.0]_, [MadDM2.0]_, [MadDM3.0]_, [MadDM3.2]_
- DRAKE [DRAKE1.0]_
- MicrOMEGAs [MicrOMEGAs2.0]_, [MicrOMEGAs3.0]_, [MicrOMEGAs4.0]_, [MicrOMEGAs5.0]_, [MicrOMEGAs5.2]_, [MicrOMEGAs5.3]_.

.. [MadDM1.0] \ M. Backovic, K. Kong, and M. McCaskey, Physics of the Dark Universe 5-6, 18 (2014), 1308.4955
.. [MadDM2.0] \ M. Backovic, A. Martini, O. Mattelaer, K. Kong, and G. Mohlabeng, Phys. Dark Univ. 9-10, 37 (2015), 1505.04190
.. [MadDM3.0] \ F. Ambrogi, C. Arina, M. Backovic, J. Heisig, F. Maltoni, L. Mantani, O. Mattelaer, and G. Mohlabeng, Phys. Dark Univ. 24, 100249 (2019), 1804.00044
.. [MadDM3.2] \ C. Arina, J. Heisig, F. Maltoni, D. Massaro, and O. Mattelaer, Eur. Phys. J. C 83, 241 (2023), 2107.04598
.. [DRAKE1.0] \ T. Binder, T. Bringmann, M. Gustafsson, and A. Hryczuk, Eur. Phys. J. C 81, 577 (2021), 2103.01944
.. [MicrOMEGAs2.0] \ G. Belanger, F. Boudjema, A. Pukhov, and A. Semenov, Comput. Phys. Commun. 176, 367 (2007), hep-ph/0607059
.. [MicrOMEGAs3.0] \ G. Belanger, F. Boudjema, A. Pukhov, and A. Semenov, Comput. Phys. Commun. 185, 960 (2014), 1305.0237
.. [MicrOMEGAs4.0] \ D. Barducci, G. Belanger, J. Bernon, F. Boudjema, J. Da Silva, S. Kraml, U. Laa, and A. Pukhov, Comput. Phys. Commun. 222, 327 (2018), 1606.03834
.. [MicrOMEGAs5.0] \ G. Belanger, F. Boudjema, A. Goudelis, A. Pukhov, and B. Zaldivar, Comput. Phys. Commun. 231, 173 (2018), 1801.03509
.. [MicrOMEGAs5.2] \ G. Belanger, A. Mjallal and A. Pukhov, Eur. Phys. J. C 81, 239 (2021), 2003.08621
.. [MicrOMEGAs5.3] \ G. Alguero, G. Belanger, S. Kraml and A. Pukhov, SciPost Phys. 13, 124 (2022), 2207.10536
