# The ABACUS Toolchain
Version 2023.2

## Author
QuantumMisaka (Zhaoqing Liu) @PKU @AISI

Inspired by cp2k-toolchain, still in improvement.

## Introduction

This toolchain will help you easily compile and install, 
or link libraries ABACUS depends on 
in ONLINE or OFFLINE way,
and give setup files that you can use to compile ABACUS.

## Todo
- [x] `gnu-openblas` toolchain support for `openmpi` and `mpich`.
- [x] `intel-mkl-mpi` toolchain support (for old oneapi using `icc` instead of `icx`).
- [x] `intelmkl-openmpi/mpich` toolchain support (need more test).
- [x] Automatic installation of [CEREAL](https://github.com/USCiLab/cereal) and [LIBNPY](https://github.com/llohse/libnpy) (by github.com)
- [ ] A better mirror station for all packages, especially for CEREAL and LIBNPY.
- [ ] A better README and Detail markdown file.
- [ ] Automatic installation of [DEEPMD](https://github.com/deepmodeling/deepmd-kit).
- [ ] Better compliation method for ABACUS-DEEPMD and ABACUS-DEEPKS.
- [ ] A better `setup` and toolchain code structure.
- [ ] Modulefile generation scripts.
- [ ] Support for [LIBRI](https://github.com/abacusmodeling/LibRI).
- [ ] Support for GPU compilation.


## Usage Online & Offline
Main script is `install_abacus_toolchain.sh`, 
which will use scripts in `scripts` directory 
to compile install dependencies of ABACUS.

```shell
> ./install_ABACUS_toolchain.sh
```

All packages will be downloaded from [cp2k-static/download](https://www.cp2k.org/static/downloads). by  `wget` , and will be detailedly compiled and installed in `install` directory by toolchain scripts, despite of `cereal` which will be downloaded from [CEREAL](https://github.com/USCiLab/cereal) and `libnpy` which will be downloaded from [LIBNPY](https://github.com/abacusmodeling/LibRI)

If one want to install ABACUS by toolchain OFFLINE, 
one can manually download all the packages and put them in `build` directory, 
then run this toolchain. 
All package will be detected and installed automatically. 
Also, one can install parts of packages OFFLINE and parts of packages ONLINE
just by using this toolchain

```shell
# for OFFLINE installation
# in toolchain directory
> mkdir build 
> cp ***.tar.gz build/
```

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

The toolchain installation process can be interrupted at anytime.
just re-run `install_abacus_toolchain.sh`, toolchain itself may fix it

If compliation is successful, a message will be shown like this:

```shell
> Done!
> To use the installed tools and libraries and ABACUS version
> compiled with it you will first need to execute at the prompt:
>   source ./install/setup
> To build ABACUS by gnu-toolchain, just use:
>     ./build_abacus_gnu.sh
> To build ABACUS by intel-toolchain, just use:
>     ./build_abacus_intel.sh
> or you can modify the builder scripts to suit your needs.
```

Then, after `source path/to/install/setup`, one can easily 
run builder scripts to build ABACUS binary software.

If users want to use toolchain but lack of some system library
dependencies, `install_requirements.sh` scripts will help.

If users want to re-install all the package, just do:
```shell
> rm -rf install/*
```
or you can also do it in a more completely way:
```shell
> rm -rf install/* build/*/* build/OpenBALS*/
```

Users can get help messages by simply:
```shell
> ./install_abacus_toolchain.sh -h # or --help
```


## More
More infomation can be read from `Details.md`, 
which is merely easily refined from cp2k-toolchain README.