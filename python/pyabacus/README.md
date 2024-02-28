Build Example: TwoCenterIntegral Section in ABACUS
==============

An example project built with [pybind11](https://github.com/pybind/pybind11)
and scikit-build-core. Python 3.7+ (see older commits for older versions of
Python).


Installation
------------

- install pybind11 and scikit-build-core by `pip install pybind11 scikit-build-core`
- clone this repository
- `pip install -v .`


CI Examples
-----------

There are examples for CI in `.github/workflows`. A simple way to produces
binary "wheels" for all platforms is illustrated in the "wheels.yml" file,
using [`cibuildwheel`][].

Use `pytest` to run all the unit tests for pyabacus in the local machine.

```shell
# pytest -v
====================================================== test session starts =======================================================
platform linux -- Python 3.8.18, pytest-8.0.0, pluggy-1.4.0 -- /root/miniconda3-gnu/envs/pyabacus/bin/python
cachedir: .pytest_cache
rootdir: /root/abacus-python/abacus-develop/python/pyabacus
configfile: pyproject.toml
testpaths: tests
collected 3 items                                                                                                                

tests/test_base_math.py::test_sphbes PASSED                                                                                [ 33%]
tests/test_base_math.py::test_sbt PASSED                                                                                   [ 66%]
tests/test_base_math.py::test_simpson PASSED                                                                               [100%]

======================================================= 3 passed in 0.14s ========================================================
```

License
-------

pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

Test call
---------

```python
import pyabacus as m
s = m.ModuleBase.Sphbes()
s.sphbesj(1, 0.0)
0.0
```

[`cibuildwheel`]:          https://cibuildwheel.readthedocs.io