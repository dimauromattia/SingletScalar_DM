# SingletScalar_DM

<!-- Put also the banners for the version of the package, the python version, the arxiv and the licence. -->

This repository contains the code for calculating cross sections, relic density
and source spectra for the singlet scalar model of dark matter.

## Installation

The package supports Python 3.10 and beyond.
Moreover, it depends on the following packages:
* [`numpy`](https://numpy.org/)
* [`scipy`](https://scipy.org/)
* [`matplotlib`](https://matplotlib.org/) (to run the examples)

To install the package along its dependencies, download or clone this repository
and run the following command in the package directory:

```shell
python3.10 -m pip install .
```

### Usage

If you have installed the package correctly, you can use it in the following
way:

```python
import singletscalar_dm as shp
```

Don't forget to [cite us!](#citation).

### Building the documentation

If, for any reason, you want to rebuild the documentation for yourself, you
need to use [`sphinx`](https://www.sphinx-doc.org/en/master/) with the configuration
present in the `docs/` folder.

You can install the required dependencies using this command in the package
root directory:

```shell
python -m pip install '.[docs]'
```

Build as usual following the [`sphinx`](https://www.sphinx-doc.org/en/master/) directives.

## Examples

The directory `examples/` contains ready-to-run examples on how to use the package
to create different plots.

Examples can be browsed as well in the [examples gallery on the github pages website](missing),
where you can also download both a python script and a jupyter notebook for each
one of them.

## Citation

If you are using this tool or the results obtained with it in a published work,
please cite [arXiv:####.#####](missing).
