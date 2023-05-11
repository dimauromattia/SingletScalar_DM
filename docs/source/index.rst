.. SingletScalar-DM documentation master file, created by
   sphinx-quickstart on Thu May 11 00:07:44 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SingletScalar-DM's documentation!
============================================

Installation
============

The package supports Python 3.10 and beyond.
Moreover, it depends on the following packages:
* .. _`numpy`: https://numpy.org/
* .. _`scipy`: https://scipy.org/
* .. _`matplotlib`: https://matplotlib.org/

To install the package along its dependencies, download a release (for stable
version) or clone this repository (for development version) and run the following
command in the package directory:

.. code-block:: bash

    python -m pip install .

Building the documentation
--------------------------

If you want to rebuild the documentation for yourself, you
need to use .. _`sphinx`: https://www.sphinx-doc.org/en/master/ with the configuration
present in the `docs/` folder.

You can install the required dependencies using this command in the package
root directory:

.. code-block:: bash

    python -m pip install '.[docs]'

Build as usual following the [`sphinx`](https://www.sphinx-doc.org/en/master/) directives.

Contents
========

.. toctree::
   Back to Github <https://github.com/dimauromattia/SingletScalar_DM>

.. toctree::
   :maxdepth: 1
   :caption: Contents:

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
