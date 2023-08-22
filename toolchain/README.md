# The ABACUS Toolchain
Version 2023.1

## Author
QuantumMisaka (Zhaoqing Liu) @PKU @AISI

Inspired by cp2k-toolchain, still in improvement.

## Introduction

This toolchain will help you easily compile and install, 
or link libraries ABACUS depends on 
and give setup files that you can use to compile ABACUS.

## Todo
- [x] `gnu-openblas` toolchain support for `openmpi` and `mpich`.
- [x] `intel-mkl-mpi` toolchain support (need more test).
- [x] `intelmkl-openmpi/mpich` toolchain support (need more test).
- [ ] Automatic installation of [CEREAL](https://github.com/USCiLab/cereal) and [LIBNPY](https://github.com/llohse/libnpy) (which is static-prepared in `lib` directory for now).
- [ ] Automatic installation of [DEEPMD](https://github.com/deepmodeling/deepmd-kit).
- [ ] Better compliation method for ABACUS-DEEPMD and ABACUS-DEEPKS.
- [ ] A better `setup` and toolchain code structure.
- [ ] Modulefile generation scripts.
- [ ] Support for [LIBRI](https://github.com/abacusmodeling/LibRI).
- [ ] Support for GPU compilation.


## Usage
Main script is `install_abacus_toolchain.sh`, 
which will use scripts in `scripts` directory 
to compile install dependencies of ABACUS.

```shell
> ./install_ABACUS_toolchain.sh
```

All packages will be downloaded from [cp2k-static/download](https://www.cp2k.org/static/downloads). by  `wget` , and will be detailedly compiled and installed in `install` directory by toolchain scripts

There are also well-modified script to run `install_abacus_toolchain.sh` for `gnu-openblas` and `intel-mkl` toolchains dependencies.

```shell
# for gnu-openblas
> ./toolchain_gnu.sh
# for intel-mkl
> ./toolchain_intel.sh
```

Users can easily compile and install dependencies of ABACUS
by running these scripts after loading `gcc` or `intel-mkl-mpi`
environment. 

If compliation is successful, a message will be shown like this:

```shell
Done!
To use the installed tools and libraries and ABACUS version
compiled with it you will first need to execute at the prompt:
  source ./install/setup
To build ABACUS by gnu-toolchain, just use:
    ./build_abacus_gnu.sh
To build ABACUS by intel-toolchain, just use:
    ./build_abacus_intel.sh
or you can modify the builder scripts to suit your needs.
```

Then, after `source path/to/install/setup`, one can easily 
run builder scripts to build ABACUS binary software.

If users want to use toolchain but lack of some system library
dependencies, `install_requirements.sh` scripts will help.

If users want to re-install all the package, just do:
```shell
rm -rf build/*/* install/* build/OpenBLAS-0.3.23/
```
or more easily:
```shell
re -rf build/ install/
```

## More
More infomation can be read from `Details.md`, 
which is merely easily refine from cp2k-toolchain README.