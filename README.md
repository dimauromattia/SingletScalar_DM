# SingletScalar-DM

[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=2305.XXXXX&color=red&style=for-the-badge)](https://arxiv.org/abs/2305.XXXXX)
[![License MIT](https://img.shields.io/static/v1?label=License&message=MIT&color=yellow&style=for-the-badge)](https://opensource.org/licenses/MIT)

This repository contains the code for calculating cross sections, relic density
and source spectra for the singlet scalar model of dark matter.

## Installation

The package supports Python 3.10 and beyond.
Moreover, it depends on the following packages:

* [`numpy`](https://numpy.org/)
* [`scipy`](https://scipy.org/)

To install the package along its dependencies, download a release (for stable
version) or clone this repository (for development version) and run the following
command in the package directory:
and run the following command in the package directory:

```shell
python3.10 -m pip install .
```
### Running the examples

The examples rely on the [`matplotlib`](https://matplotlib.org/) package in order to show the plots.
You can install this dependence using the following command in the package root
directory:

```shell
python -m pip install '.[examples]'
```

### Building the documentation

If you want to rebuild the documentation for yourself, you need to use [`sphinx`](https://www.sphinx-doc.org/en/master/)
with the configuration provided in the `docs/` folder.

You can install the required dependencies using the following command in the package
root directory:

```shell
python -m pip install '.[docs]'
```

Build as usual following the [`sphinx`](https://www.sphinx-doc.org/en/master/).

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
please cite [arXiv:2305.#####](https://arxiv.org/abs/2305.XXXXX).
