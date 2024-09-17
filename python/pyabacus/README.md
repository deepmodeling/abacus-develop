Build Example: TwoCenterIntegral Section in ABACUS
==================================================

An example project built with [pybind11](https://github.com/pybind/pybind11)
and scikit-build-core. Python 3.7+ (see older commits for older versions of
Python).

Installation
------------

- Create and activate a new conda env, e.g. `conda create -n myenv python=3.8 & conda activate myenv`.
- Clone ABACUS main repository and `cd abacus-develop/python/pyabacus`.
- Build pyabacus by `pip install -v .` or install test dependencies & build  pyabacus by `pip install .[test]`. (Use `pip install -v .[test] -i https://pypi.tuna.tsinghua.edu.cn/simple` to accelerate installation process.)

CI Examples
-----------

There are examples for CI in `.github/workflows`. A simple way to produces
binary "wheels" for all platforms is illustrated in the "wheels.yml" file,
using .

Use `pytest -v` to run all the unit tests for pyabacus in the local machine.

```shell
$ cd tests/
$ pytest -v
```

Run `python vis_nao.py` to visualize the numerical orbital.

```shell
$ cd examples/
$ python vis_nao.py
```

Run `python ex_s_rotate.py` in `examples` to check the S matrix.

```shell
$ cd examples/
$ python ex_s_rotate.py
norm(S_e3 - S_numer) =  3.341208104032616e-15
```

Run `python diago_matrix.py` in `examples` to check the diagonalization of a matrix.

```shell
$ cd examples/
$ python diago_matrix.py
====== Calculating eigenvalues using davidson method... ======
eigenvalues calculated by pyabacus-davidson is: 
 [-0.38440611  0.24221155  0.31593272  0.53144616  0.85155108  1.06950154
  1.11142051  1.1246216 ]
eigenvalues calculated by scipy:  [-0.38440611  0.24221155  0.31593272  0.53144616  0.85155108  1.06950154
  1.11142051  1.12462151]
eigenvalues difference:  [4.25992575e-13 2.58706945e-12 8.26848034e-11 9.64006652e-13
 1.38805634e-11 1.61699987e-10 1.05329057e-09 8.97058461e-08]

====== Calculating eigenvalues using dav_subspace method... ======
enter diag... is_subspace = 0, ntry = 0
eigenvalues calculated by pyabacus-dav_subspace is: 
 [-0.38440611  0.24221155  0.31593272  0.53144616  0.85155108  1.06950154
  1.11142052  1.12462151]
eigenvalues calculated by scipy:  [-0.38440611  0.24221155  0.31593272  0.53144616  0.85155108  1.06950154
  1.11142051  1.12462151]
eigenvalues difference:  [4.98749930e-11 5.52219381e-12 1.05679354e-11 3.20832250e-12
 4.96347408e-12 7.22339299e-10 4.29339986e-09 6.92761959e-09]
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

[`cibuildwheel`]: https://cibuildwheel.readthedocs.io
