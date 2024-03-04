Build Example: TwoCenterIntegral Section in ABACUS
==============

An example project built with [pybind11](https://github.com/pybind/pybind11)
and scikit-build-core. Python 3.7+ (see older commits for older versions of
Python).


Installation
------------

- Create and activate a new conda env, e.g. `conda create -n myenv python=3.8 & conda activate myenv`.
- Clone ABACUS main repository and `cd abacus-develop/python/pyabacus`.
- Install dependencies by `pip install -r requirements.txt`.(Use `pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple` to accelerate installation process.)
- Build pyabacus by `pip install -v .`.


CI Examples
-----------

There are examples for CI in `.github/workflows`. A simple way to produces
binary "wheels" for all platforms is illustrated in the "wheels.yml" file,
using [`cibuildwheel`][].

Use `pytest` to run all the unit tests for pyabacus in the local machine.

Run `python test_m_nao.py` in `tests` to check the S matrix.

```shell
$ python test_m_nao.py 
norm(S_e3 - S_numer) =  3.341208104032616e-15
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