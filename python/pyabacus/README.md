Build Example: TwoCenterIntegral Section in ABACUS
==============

An example project built with [pybind11](https://github.com/pybind/pybind11)
and scikit-build-core. Python 3.7+ (see older commits for older versions of
Python).


Installation
------------

- install pybind11 and scikit-build-core by `pip install pybind11 scikit-build-core`
- clone this repository
- `pip install ./pyabacus`


CI Examples
-----------

There are examples for CI in `.github/workflows`. A simple way to produces
binary "wheels" for all platforms is illustrated in the "wheels.yml" file,
using [`cibuildwheel`][].

License
-------

pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

Test call
---------

```python
import pyabacus as m

chi = m.NumericalRadial()
```

[`cibuildwheel`]:          https://cibuildwheel.readthedocs.io