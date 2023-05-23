# SingletScalar-DM

[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=2305.11937&color=red&style=for-the-badge)](https://arxiv.org/abs/2305.11937)
[![License MIT](https://img.shields.io/static/v1?label=License&message=MIT&color=yellow&style=for-the-badge)](https://opensource.org/licenses/MIT)

This repository contains the code for calculating cross sections, relic density
and source spectra for the singlet scalar model of dark matter.

The package relies on pre-computed data obtained according to several codes:
MadDM \[[1](#maddm-1)&ndash;[4](#maddm-4)\], DRAKE \[[5](#drake)\] and micrOMEGAs \[[6](#micromegas-1)&ndash;[11](#micromegas-6)\].

## Installation

The package supports Python 3.7 and beyond.
Moreover, it depends on the following packages:

* [`numpy`](https://numpy.org/)
* [`scipy`](https://scipy.org/)

To install the package along its dependencies, download a release (for stable
version) or clone this repository (for development version) and run the following
command in the package directory:
and run the following command in the package directory:

```bash
python -m pip install .
```
### Running the examples

The examples rely on the [`matplotlib`](https://matplotlib.org/) package in order to show the plots.
You can install this dependence using the following command in the package root
directory:

```bash
python -m pip install '.[examples]'
```

### Building the documentation

If you want to rebuild the documentation for yourself, you need to use [`sphinx`](https://www.sphinx-doc.org/en/master/)
with the configuration provided in the `docs/` folder.

You can install the required dependencies using the following command in the package
root directory:

```bash
python -m pip install '.[docs]'
```

Build using [`sphinx`](https://www.sphinx-doc.org/en/master/).
You can run the following commands inside the `docs/` directory.

If you want to update the source documentation:

```bash
make rst
```

If you want to build the documentation in a certain format, you can specify a 
``<builder>`` among the [available builders](https://www.sphinx-doc.org/en/master/usage/builders/index.html):

```bash
make <builder>
```

Notice that it can take a while to build the whole documentation, given it will
run all the examples.

### Install in development mode

If you wish to change something in the data files or in the code, we advise you
to install the package in *development* or *editable* mode, using the following
command in the package root directory:

```bash
python -m pip install --editable .
```

In this way you are free to modify anything inside the `src/` folder and the changes
will be immediately reflected in the installed package once you re-import it.

### Usage

If you have installed the package correctly, you can use it in the following
way:

```python
from singletscalar_dm import *
```

Don't forget to [cite us!](#citation).

## Examples

The directory `examples/` contains ready-to-run examples on how to use the package
to create different plots.
Refer to [this section](#running-the-examples) to install the required dependencies.

Examples can be browsed as well in the [examples gallery on the GitHub pages website](missing),
where you can also download both a python script and a jupyter notebook for each
one of them.

## Citation

If you are using this tool or the results obtained with it in a published work,
please cite [arXiv:2305.11937](https://arxiv.org/abs/2305.11937).

## References

\[<a id="maddm-1">1</a>\] M. Backovic, K. Kong, and M. McCaskey, Physics of the Dark Universe 5-6, 18 (2014), 1308.4955  
\[2\] M. Backovic, A. Martini, O. Mattelaer, K. Kong, and G. Mohlabeng, Phys. Dark Univ. 9-10, 37 (2015), 1505.04190  
\[3\] F. Ambrogi, C. Arina, M. Backovic, J. Heisig, F. Maltoni, L. Mantani, O. Mattelaer, and G. Mohlabeng, Phys. Dark Univ. 24, 100249 (2019), 1804.00044  
\[<a id="maddm-4">4</a>\] C. Arina, J. Heisig, F. Maltoni, D. Massaro, and O. Mattelaer, Eur. Phys. J. C 83, 241 (2023), 2107.04598  
\[<a id="drake">5</a>\] T. Binder, T. Bringmann, M. Gustafsson, and A. Hryczuk, Eur. Phys. J. C 81, 577 (2021), 2103.01944  
\[<a id="micromegas-1">6</a>\] G. Belanger, F. Boudjema, A. Pukhov, and A. Semenov, Comput. Phys. Commun. 176, 367 (2007), hep-ph/0607059  
\[7\] G. Belanger, F. Boudjema, A. Pukhov, and A. Semenov, Comput. Phys. Commun. 185, 960 (2014), 1305.0237  
\[8\] D. Barducci, G. Belanger, J. Bernon, F. Boudjema, J. Da Silva, S. Kraml, U. Laa, and A. Pukhov, Comput. Phys. Commun. 222, 327 (2018), 1606.03834  
\[9\] G. Belanger, F. Boudjema, A. Goudelis, A. Pukhov, and B. Zaldivar, Comput. Phys. Commun. 231, 173 (2018), 1801.03509  
\[10\] G. Belanger, A. Mjallal and A. Pukhov, Eur. Phys. J. C 81, 239 (2021), 2003.08621  
\[<a id="micromegas-6">11</a>\] G. Alguero, G. Belanger, S. Kraml and A. Pukhov, SciPost Phys. 13, 124 (2022), 2207.10536  
