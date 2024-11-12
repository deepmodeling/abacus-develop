# The ABACUS Toolchain

Version 2024.3

## Author

[QuantumMisaka](https://github.com/QuantumMisaka) 
(Zhaoqing Liu) @PKU @AISI

Inspired by cp2k-toolchain, still in improvement.

You should have read this README before using this toolchain.

## Introduction

This toolchain will help you easily compile and install, 
or link libraries ABACUS depends on 
in ONLINE or OFFLINE way,
and give setup files that you can use to compile ABACUS.

## Todo

- [x] `gnu-openblas` toolchain support for `openmpi` and `mpich`.
- [x] `intel-mkl-mpi` toolchain support using `icc`/`icpc`/`ifort` or `icx`/`icpx`/`ifort`. (`icx` as default, but will have problem for ELPA in AMD machine, one can specify `--with-intel-classic=yes` to use `icc`), 
- [x] `intel-mkl-mpich` toolchain support.
- [x] Automatic installation of [CEREAL](https://github.com/USCiLab/cereal) and [LIBNPY](https://github.com/llohse/libnpy) (by github.com)
- [x] Support for [LibRI](https://github.com/abacusmodeling/LibRI) by submodule or automatic installation from github.com (but installed LibRI via `wget` seems to have some problem, please be cautious)
- [x] A mirror station by Bohrium database, which can download CEREAL, LibNPY, LibRI and LibComm by `wget` in China Internet. 
- [x] Support for GPU compilation, users can add `-DUSE_CUDA=1` in builder scripts.
- [ ] Change the downloading url from cp2k mirror to other mirror or directly downloading from official website. (doing)
- [ ] A better README and Detail markdown file.
- [ ] Automatic installation of [DEEPMD](https://github.com/deepmodeling/deepmd-kit).
- [ ] Better compliation method for ABACUS-DEEPMD and ABACUS-DEEPKS.
- [ ] Modulefile generation scripts.
- [ ] Support for AMD compiler and math lib like `AOCL` and `AOCC`


## Usage Online & Offline

Main script is *install_abacus_toolchain.sh*, 
which will use scripts in *scripts* directory 
to compile install dependencies of ABACUS.
It can be directly used, but not recommended.

There are also well-modified script to run *install_abacus_toolchain.sh* for `gnu-openblas` and `intel-mkl` toolchains dependencies.

```shell
# for gnu-openblas
> ./toolchain_gnu.sh
# for intel-mkl
> ./toolchain_intel.sh
# for intel-mkl-mpich
> ./toolchain_intel-mpich.sh
```

It is recommended to run one of them first to get a fast installation of ABACUS under certain environments.

If you are using Intel environments via Intel-OneAPI: please note:
1. After version 2024.0, Intel classic compilers `icc` and `icpc` are not present, so as `ifort` after version 2025.0. Intel MPI compiler will also be updated to `mpiicx`, `mpiicpx` and `mpiifx`.
2. toolchain will detect `icx`, `icpx`, `ifx`, `mpiicx`, `mpiicpx` and `mpiifx` as default compiler.
3. Users can manually specify `--with-intel-classic=yes` to use Intel classic compiler in `toolchain*.sh`, or specify `--with-intel-mpi-clas=yes` to use Intel MPI classic compiler in `toolchain*.sh` while keep the CC, CXX and F90 compiler to new version.
4. Users can manually specify `--with-ifx=no` in `toolchain*.sh` to use `ifort` while keep other compiler to new version. 
5. More information is in the later part of this README.

**Notice: You GCC version should be no lower than 5 !!!, larger than 7.3.0 is recommended**

**Notice: You SHOULD `source` or `module load` related environments before use toolchain method for installation, espacially for `gcc` or `intel-oneAPI` !!!! for example, `module load mkl mpi icc compiler`**

**Notice: You SHOULD keep your environments systematic, for example, you CANNOT load `intel-OneAPI` environments while use gcc toolchain !!!**

**Notice: If your server system already have libraries like `cmake`, `openmpi`, please change related setting in `toolchain*.sh` like `--with-cmake=system`**


All packages will be downloaded from [cp2k-static/download](https://www.cp2k.org/static/downloads). by  `wget` , and will be detailedly compiled and installed in `install` directory by toolchain scripts, despite of:

- `CEREAL` which will be downloaded from [CEREAL](https://github.com/USCiLab/cereal)  
- `Libnpy` which will be downloaded from [LIBNPY](https://github.com/llohse/libnpy)
- `LibRI` which will be downloaded from [LibRI](https://github.com/abacusmodeling/LibRI)
- `LibCOMM` which will be downloaded from [LibComm](https://github.com/abacusmodeling/LibComm)
- `RapidJSON` which will be downloaded from [RapidJSON](https://github.com/Tencent/rapidjson)
Notice: These packages will be downloaded by `wget` from `github.com`, which is hard to be done in Chinese Internet. You may need to use offline installation method. 

Instead of github.com, we offer other package station, you can use it by:
```shell
wget https://bohrium-api.dp.tech/ds-dl/abacus-deps-93wi-v3 -O abacus-deps-v3.zip
```
`unzip` it ,and you can do offline installation of these packages above after rename. 
```shell
# packages downloaded from github.com
mv v1.3.2.tar.gz build/cereal-1.3.2.tar.gz
```
The above station will be updated handly but one should notice that the version will always lower than github repo.

If one want to install ABACUS by toolchain OFFLINE, 
one can manually download all the packages from [cp2k-static/download](https://www.cp2k.org/static/downloads) or official website
and put them in *build* directory by formatted name
like *fftw-3.3.10.tar.gz*, or *openmpi-5.0.5.tar.bz2*, 
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

The needed dependencies version default:

- `cmake` 3.30.0
- `gcc` 13.2.0 (which will always NOT be installed, But use system)
- `OpenMPI` 4.1.6 (5.0.5 can be used but have some problem in OpenMP parallel computation in ELPA)
- `MPICH` 4.2.2
- `OpenBLAS` 0.3.28 (Intel toolchain need `get_vars.sh` tool from it)
- `ScaLAPACK` 2.2.1 (a developing version)
- `FFTW` 3.3.10
- `LibXC` 6.2.2
- `ELPA` 2024.05.001
- `CEREAL` 1.3.2
- `RapidJSON` 1.1.0
And Intel-oneAPI need user or server manager to manually install from Intel.
[Intel-oneAPI](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/toolkits.html)

Dependencies below are optional， which is NOT installed by default:

- `LibTorch` 2.1.2
- `Libnpy` 1.0.1
- `LibRI` 0.2.0
- `LibComm` 0.1.1

Users can install them by using `--with-*=install` in toolchain*.sh, which is `no` in default.
> Notice: LibRI, LibComm and Libnpy is on actively development, you should check-out the package version when using this toolchain. Also, LibRI and LibComm can be installed by github submodule, that is also work for libnpy, which is more recommended.

Users can easily compile and install dependencies of ABACUS
by running these scripts after loading `gcc` or `intel-mkl-mpi`
environment. 

The toolchain installation process can be interrupted at anytime.
just re-run *toolchain_\*.sh*, toolchain itself may fix it

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

You can run *build_abacus_gnu.sh* or *build_abacus_intel.sh* to build ABACUS 
by gnu-toolchain or intel-toolchain respectively, the builder scripts will
automatically locate the environment and compile ABACUS.
You can manually change the builder scripts to suit your needs.
The builder scripts will generate `abacus_env.sh` for source

Then, after `source abacus_env.sh`, one can easily 
run builder scripts to build ABACUS binary software.

If users want to use toolchain but lack of some system library
dependencies, *install_requirements.sh* scripts will help.

If users want to re-install all the package, just do:

```shell
> rm -rf install
```

or you can also do it in a more completely way:

```shell
> rm -rf install build/*/* build/OpenBLAS*/ build/setup_*
```

## Common Problems and Solutions

### LibRI and LibComm for EXX

- GCC toolchain with OpenMPI cannot compile LibComm v0.1.1 due to the different MPI variable type from MPICH and IntelMPI, see discussion here [#5033](https://github.com/deepmodeling/abacus-develop/issues/5033), you can switch to GCC-MPICH or Intel toolchain
- It is recommended to use Intel toolchain if one wants to include EXX feature in ABACUS, which can have much better performance and can use more than 16 threads in OpenMP parallelization to accelerate the EXX process.

### GPU version of ABACUS

For GPU version of ABACUS (do not GPU version installer of ELPA, which is still doing work), add following options in build*.sh:

```shell
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=icpx \
        -DMPI_CXX_COMPILER=mpiicpc \
        ......
        -DUSE_CUDA=1 \
        -DCMAKE_CUDA_COMPILER=${path to cuda toolkit}/bin/nvcc \
        ......
```

Notice: You CANNOT use `icpx` compiler for GPU version of ABACUS for now, see discussion here [#2906](https://github.com/deepmodeling/abacus-develop/issues/2906) and [#4976](https://github.com/deepmodeling/abacus-develop/issues/4976)

If you wants to use ABACUS GPU-LCAO by `cusolvermp` or `elpa`, please contact the coresponding developer, toolchain do not fully support them now.

### Shell problem

If you encounter problem like:

```shell
/bin/bash^M: bad interpreter: No such file or directory
```

or   `permission denied` problem, you can simply run:

```shell
./pre_set.sh
```

And also, you can fix `permission denied` problem via `chmod +x`
if *pre_set.sh* have no execution permission; 
if the *pre_set.sh* also have `/bin/bash^M` problem, you can run:

```shell
> dos2unix pre_set.sh
```

to fix it

### Libtorch and DeePKS problem

If deepks feature have problem, you can manually change libtorch version
from 2.1.2 to 2.0.1 or 1.12.0 in `toolchain/scripts/stage4/install_libtorch.sh`.

Also, you can install ABACUS without deepks by removing all the deepks and related options.

NOTICE: if you want deepks feature, your intel-mkl environment should be accessible in building process. you can check it in `build_abacus_gnu.sh`

### DeePMD feature problem

When you encounter problem like `GLIBCXX_3.4.29 not found`, it is sure that your `gcc` version is lower than the requirement of `libdeepmd`.

After my test, you need `gcc`>11.3.1 to enable deepmd feature in ABACUS.

### Intel-oneAPI problem

#### ELPA problem via Intel-oneAPI toolchain in AMD server

The default compiler for Intel-oneAPI is `icpx` and `icx`, which will cause problem when compling ELPA in AMD server. (Which is a problem and needed to have more check-out)

The best way is to change `icpx` to `icpc`, `icx` to `icc`. user can manually change it in toolchain*.sh via `--with-intel-classic=yes`

Notice: `icc` and `icpc` from Intel Classic Compiler of Intel-oneAPI is not supported for 2024.0 and newer version. And Intel-OneAPI 2023.2.0 can be found in website. See discussion here [#4976](https://github.com/deepmodeling/abacus-develop/issues/4976)

#### link problem in early 2023 version oneAPI

Sometimes Intel-oneAPI have problem to link `mpirun`, 
which will always show in 2023.2.0 version of MPI in Intel-oneAPI. 
Try `source /path/to/setvars.sh` or install another version of IntelMPI may help.

which is fixed in 2024.0.0 version of Intel-oneAPI, 
And will not occur in Intel-MPI before 2021.10.0 (Intel-oneAPI before 2023.2.0)

More problem and possible solution can be accessed via [#2928](https://github.com/deepmodeling/abacus-develop/issues/2928)

## Advanced Installation Usage

1. Users can move toolchain directory to anywhere you like, 
and complete installation by change the setting in 
`toolchain_*.sh` and `build_*.sh` by your own setting.
By moving toolchain out or rename it ,one can make toolchain independent
from ABACUS repo, make dependencies package more independent and flexible.
2. Users can manually change `pkg_install_dir` variable 
in `scripts/stage*/install*` to change the installation directory 
of each packages, which may let the installation more fiexible.


## More

More infomation can be read from `Details.md`.