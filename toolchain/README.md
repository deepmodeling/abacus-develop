# The ABACUS Toolchain

Version 2025.1

## Main Developer

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
- [x] Support for GPU-PW and GPU-LCAO compilation (elpa, cusolvermp is developing), and `-DUSE_CUDA=1` is needed builder scripts.
- [x] Support for AMD compiler and math lib  `AOCL` and `AOCC` (not fully complete due to flang and AOCC-ABACUS compliation error)
- [ ] Support for more GPU device out of Nvidia.
- [ ] Change the downloading url from cp2k mirror to other mirror or directly downloading from official website. (doing)
- [ ] Support a JSON or YAML configuration file for toolchain, which can be easily modified by users.
- [ ] A better README and Detail markdown file.
- [ ] Automatic installation of [DEEPMD](https://github.com/deepmodeling/deepmd-kit).
- [ ] Better compliation method for ABACUS-DEEPMD and ABACUS-DEEPKS.
- [ ] Modulefile generation scripts.


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
# for amd aocc-aocl
> ./toolchain_amd.sh
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
like *fftw-3.3.10.tar.gz*, or *openmpi-5.0.6.tar.bz2*, 
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

- `cmake` 3.31.2
- `gcc` 13.2.0 (which will always NOT be installed, But use system)
- `OpenMPI` 5.0.6 (Version 5 OpenMPI is good but will have compability problem, user can manually downarade to Version 4 in toolchain scripts)
- `MPICH` 4.3.0
- `OpenBLAS` 0.3.28 (Intel toolchain need `get_vars.sh` tool from it)
- `ScaLAPACK` 2.2.1 (a developing version)
- `FFTW` 3.3.10
- `LibXC` 7.0.0
- `ELPA` 2025.01.001
- `CEREAL` 1.3.2
- `RapidJSON` 1.1.0
And:
- Intel-oneAPI need user or server manager to manually install from Intel.
- - [Intel-oneAPI](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/toolkits.html)
- AMD AOCC-AOCL need user or server manager to manually install from AMD.
- - [AOCC](https://www.amd.com/zh-cn/developer/aocc.html)
- - [AOCL](https://www.amd.com/zh-cn/developer/aocl.html)

Dependencies below are optional， which is NOT installed by default:

- `LibTorch` 2.1.2
- `Libnpy` 1.0.1
- `LibRI` 0.2.0
- `LibComm` 0.1.1

Users can install them by using `--with-*=install` in toolchain*.sh, which is `no` in default. Also, user can specify the absolute path of the package by `--with-*=path/to/package` in toolchain*.sh to allow toolchain to use the package.
> Notice: LibTorch always suffer from GLIBC_VERSION problem, if you encounter this, please downgrade LibTorch version to 1.12.1 in scripts/stage4/install_torch.sh
> 
> Notice: LibRI, LibComm, Rapidjson and Libnpy is on actively development, you should check-out the package version when using this toolchain. 

Users can easily compile and install dependencies of ABACUS
by running these scripts after loading `gcc` or `intel-mkl-mpi`
environment. 

The toolchain installation process can be interrupted at anytime.
just re-run *toolchain_\*.sh*, toolchain itself may fix it. If you encouter some problem, you can always remove some package in the interrupted points and re-run the toolchain.

Some useful options:
- `--dry-run`: just run the main install scripts for environment setting, without any package downloading or installation.
- `--pack-run`: just run the install scripts without any package building, which helps user to download and check the packages, paticularly for offline installation to a server.

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
> To build ABACUS by amd-toolchain in gcc-aocl, just use:
>     ./build_abacus_amd.sh
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

## GPU version of ABACUS

Toolchain supports compiling GPU version of ABACUS with Nvidia-GPU and CUDA. For usage, adding following options in build*.sh:

```shell
# in build_abacus_gnu.sh
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=g++ \
        -DMPI_CXX_COMPILER=mpicxx \
        ......
        -DUSE_CUDA=ON \
        # -DCMAKE_CUDA_COMPILER=${path to cuda toolkit}/bin/nvcc \ # add if needed
        ......
# in build_abacus_intel.sh
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=icpc \
        -DMPI_CXX_COMPILER=mpiicpc \
        ......
        -DUSE_CUDA=ON \
        # -DCMAKE_CUDA_COMPILER=${path to cuda toolkit}/bin/nvcc \ # add if needed
        ......
```
which will enable GPU version of ABACUS, and the `ks_solver cusolver` method can be directly used for PW and LCAO calculation.

Notice: You CANNOT use `icpx` compiler for GPU version of ABACUS for now, see discussion here [#2906](https://github.com/deepmodeling/abacus-develop/issues/2906) and [#4976](https://github.com/deepmodeling/abacus-develop/issues/4976)

If you wants to use ABACUS GPU-LCAO by `cusolvermp` or `elpa` for multiple-GPU calculation, please compile according to the following usage:

1. For the elpa method, add
```shell
export CUDA_PATH=/path/to/CUDA
# install_abacus_toolchain.sh part options
--enable-cuda \
--gpu-ver=(GPU-compatibility-number) \
```
to the `toolchain_*.sh`, and then follow the normal step to install the dependencies using `./toolchain_*.sh`. For checking the GPU compatibility number, you can refer to the [CUDA compatibility](https://developer.nvidia.com/cuda-gpus).

Afterwards, make sure these option are enable in your `build_abacus_*.sh` script 
```shell
-DUSE_ELPA=ON \
-DUSE_CUDA=ON \
```
then just build the abacus executable program by compiling it with `./build_abacus_*.sh`.

The ELPA method need more parameter setting, but it doesn't seem to be affected by the CUDA toolkits version, and it is no need to manually install and package. 

2. For the cusolvermp method, toolchain_*.sh does not need to be changed, just follow it directly install dependencies using `./toolchain_*.sh`, and then add
```shell
-DUSE_CUDA=ON \
-DENABLE_CUSOLVERMP=ON \
-D CAL_CUSOLVERMP_PATH=/path/to/math.libs/1x.x/target/x86_64-linux/lib \
```
to the `build.abacus_*.sh` file. add the following three items to the environment (assuming you are using hpcsdk):
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/comm_libs/1x.x/hpcx/hpcx-x.xx/ucc/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/comm_libs/1x.x/hpcx/hpcx-x.xx/ucx/lib
export CPATH=$CPATH:/path/to/math_libs/1x.x/targets/x86_64-linux/include
```
Just enough to build the abacus executable program by compiling it with `./build_abacus_*.sh`.

You can refer to the linking video for auxiliary compilation and installation. [Bilibili](https://www.bilibili.com/video/BV1eqr5YuETN/).

The cusolverMP requires installation from sources such as apt or yum, which is suitable for containers or local computers.
The second choice is using [NVIDIA HPC_SDK](https://developer.nvidia.com/hpc-sdk-downloads) for installation, which is relatively simple, but the package from NVIDIA HPC_SDK may not be suitable, especially for muitiple-GPU parallel running. To better use cusolvermp and its dependency (libcal, ucx, ucc) in multi-GPU running, please contact your server manager.

After compiling, you can specify `device GPU` in INPUT file to use GPU version of ABACUS.


## Common Problems and Solutions

### Intel-oneAPI problem

#### OneAPI 2025.0 problem

Generally, OneAPI 2025.0 can be useful to compile basic function of ABACUS, but one will encounter compatible problem related to something. Here is the treatment
- related to rapidjson: 
- - Not to use rapidjson in your toolchain
- - or use the master branch of [RapidJSON](https://github.com/Tencent/rapidjson)
- related to LibRI: not to use LibRI or downgrade your OneAPI.

#### ELPA problem via Intel-oneAPI toolchain in AMD server

The default compiler for Intel-oneAPI is `icpx` and `icx`, which will cause problem when compling ELPA in AMD server. (Which is a problem and needed to have more check-out)

The best way is to change `icpx` to `icpc`, `icx` to `icc`. user can manually change it in *toolchain_intel.sh* via `--with-intel-classic=yes`

Notice: `icc` and `icpc` from Intel Classic Compiler of Intel-oneAPI is not supported for 2024.0 and newer version. And Intel-OneAPI 2023.2.0 can be found in QE website. You need to download Base-toolkit for MKL and HPC-toolkit for MPi and compiler for Intel-OneAPI 2023.2.0, while in Intel-OneAPI 2024.x, only the HPC-toolkit is needed.

You can get Intel-OneAPI in [QE-managed website](https://pranabdas.github.io/espresso/setup/hpc/#installing-intel-oneapi-libraries), and use this code to get Intel oneAPI Base Toolkit and HPC Toolkit:
```shell
wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397_offline.sh
wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/0722521a-34b5-4c41-af3f-d5d14e88248d/l_HPCKit_p_2023.2.0.49440_offline.sh
```

Related discussion here [#4976](https://github.com/deepmodeling/abacus-develop/issues/4976)

#### linking problem in early 2023 version oneAPI

Sometimes Intel-oneAPI have problem to link `mpirun`, 
which will always show in 2023.2.0 version of MPI in Intel-oneAPI. 
Try `source /path/to/setvars.sh` or install another version of IntelMPI may help.

which is fixed in 2024.0.0 version of Intel-oneAPI, 
And will not occur in Intel-MPI before 2021.10.0 (Intel-oneAPI before 2023.2.0)

More problem and possible solution can be accessed via [#2928](https://github.com/deepmodeling/abacus-develop/issues/2928)

### AMD AOCC-AOCL problem

You cannot use AOCC to complie abacus now, see [#5982](https://github.com/deepmodeling/abacus-develop/issues/5982) .

However, use AOCC-AOCL to compile dependencies is permitted and usually get boosting in ABACUS effciency. But you need to get rid of `flang` while compling ELPA. Toolchain itself help you make this `flang` shade in default, and you can manully use `flang` by setting `--with-flang=yes` in `toolchain_amd.sh` to have a try. 

Notice: ABACUS via GCC-AOCL in AOCC-AOCL toolchain have no application with DeePKS, DeePMD and LibRI. 

### OpenMPI problem

#### in EXX and LibRI

- GCC toolchain with OpenMPI cannot compile LibComm v0.1.1 due to the different MPI variable type from MPICH and IntelMPI, see discussion here [#5033](https://github.com/deepmodeling/abacus-develop/issues/5033), you can try use a newest branch of LibComm by 
```
git clone https://gitee.com/abacus_dft/LibComm -b MPI_Type_Contiguous_Pool
``` 
or pull the newest master branch of LibComm
```
git clone https://github.com/abacusmodeling/LibComm
```
. yet another is switching to GCC-MPICH or Intel toolchain
- It is recommended to use Intel toolchain if one wants to include EXX feature in ABACUS, which can have much better performance and can use more than 16 threads in OpenMP parallelization to accelerate the EXX process.

#### OpenMPI-v5 

OpenMPI in version 5 has huge update, lead to compatibility problem. If one wants to use the OpenMPI in version 4 (4.1.6), one can specify `--with-openmpi-4th=yes` in *toolchain_gnu.sh*


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
