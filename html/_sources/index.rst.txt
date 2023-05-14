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
