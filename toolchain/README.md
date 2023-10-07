# The ABACUS Toolchain
Version 2023.4

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
- [x] `intel-mkl-mpi` toolchain support using `icc` or `icx`. (`icx` version of ABACUS have some problem now, please be cautious)
- [x] `intel-mkl-mpich` toolchain support.
- [x] Automatic installation of [CEREAL](https://github.com/USCiLab/cereal) and [LIBNPY](https://github.com/llohse/libnpy) (by github.com)
- [x] Support for [LibRI](https://github.com/abacusmodeling/LibRI) by submodule or automatic installation from github.com (but installed LibRI via `wget` seems to have some problem, please be cautious)
- [x] A mirror station by Bohrium database, which can download CEREAL, LibNPY, LibRI and LibComm by `wget` in China Internet. 
- [x] Support for GPU compilation, users can add `-DUSE_CUDA=1` in builder scripts.
- [ ] A better mirror station for all packages.
- [ ] A better README and Detail markdown file.
- [ ] Automatic installation of [DEEPMD](https://github.com/deepmodeling/deepmd-kit).
- [ ] Better compliation method for ABACUS-DEEPMD and ABACUS-DEEPKS.
- [ ] A better `setup` and toolchain code structure.
- [ ] Modulefile generation scripts.
- [ ] Support for `acml` toolchain (scripts are partly in toolchain now) or other AMD compiler and math lib like `AOCL` and `AOCC`


## Usage Online & Offline
Main script is `install_abacus_toolchain.sh`, 
which will use scripts in `scripts` directory 
to compile install dependencies of ABACUS.

**Notice: You SHOULD `source` or `module load` related environments before use toolchain method for installation, espacially for `gcc` or `intel-oneAPI` !!!! for example, `module load mkl mpi icc compiler`**

**Notice: You SHOULD keep your environments systematic, for example, you CANNOT load `intel-OneAPI` environments while use gcc toolchain !!!**

**Notice: If your server system already have libraries like `cmake`, `openmpi`, please change related setting in `toolchain*.sh` like `--with-cmake=system`**

```shell
> ./install_ABACUS_toolchain.sh
```

All packages will be downloaded from [cp2k-static/download](https://www.cp2k.org/static/downloads). by  `wget` , and will be detailedly compiled and installed in `install` directory by toolchain scripts, despite of:
- `CEREAL` which will be downloaded from [CEREAL](https://github.com/USCiLab/cereal)  
- `LibNPY` which will be downloaded from [LIBNPY](https://github.com/llohse/libnpy)
- `LibRI` which will be downloaded from [LibRI](https://github.com/abacusmodeling/LibRI)
- `LibCOMM` which will be downloaded from [LibComm](https://github.com/abacusmodeling/LibComm)
Notice: These packages will be downloaded by `wget` from `github.com`, which is hard to be done in Chinese Internet. You may need to use offline installation method. 

Instead of github.com, we offer other package station, you can use it by:
```shell
wget https://bohrium-api.dp.tech/ds-dl/abacus-deps-93wi-v1 -O abacus-deps-v1.zip
```
`unzip` it ,and you can do offline installation of these packages above

If one want to install ABACUS by toolchain OFFLINE, 
one can manually download all the packages from [cp2k-static/download](https://www.cp2k.org/static/downloads) or official website
and put them in `build` directory by formatted name
like `fftw-3.3.10.tar.gz`, or `openmpi-4.1.5.tar.gz`, 
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

Notice: for `CEREAL`, `LibNPY`, `LibRI` and `LibCOMM`, 
you need to download them from github.com, 
rename it as formatted, and put them in `build` directory at the same time
e.g.:
```shell
# packages downloaded from github.com
mv v1.3.2.tar.gz build/cereal-1.3.2.tar.gz
```

There are also well-modified script to run `install_abacus_toolchain.sh` for `gnu-openblas` and `intel-mkl` toolchains dependencies.

```shell
# for gnu-openblas
> ./toolchain_gnu.sh
# for intel-mkl
> ./toolchain_intel.sh
# for intel-mkl-mpich
> ./toolchain_intel-mpich.sh
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
You can run build_abacus_gnu.sh or build_abacus_intel.sh to build ABACUS 
by gnu-toolchain or intel-toolchain respectively, the builder scripts will
automatically locate the environment and compile ABACUS.
You can manually change the builder scripts to suit your needs.
The builder scripts will generate `abacus_env.sh` for source

Then, after `source abacus_env.sh`, one can easily 
run builder scripts to build ABACUS binary software.

If users want to use toolchain but lack of some system library
dependencies, `install_requirements.sh` scripts will help.

If users want to re-install all the package, just do:
```shell
> rm -rf install
```
or you can also do it in a more completely way:
```shell
> rm -rf install build/*/* build/OpenBLAS*/ build/setup_*
```

Users can get help messages by simply:
```shell
> ./install_abacus_toolchain.sh -h # or --help
```


## Common Problem and Solution
If you encounter problem like:
```shell
/bin/bash^M: bad interpreter: No such file or directory
```
or   `permission denied` problem, you can simply run:
```shell
./pre_set.sh
```
And also, you can fix `permission denied` problem via `chmod +x`
if `pre_set.sh` have no execution permission.


## Advanced Installation Usage

1. Users can move toolchain directory to anywhere you like, 
and complete installation by change the setting in 
`toolchain_*.sh` and `build_*.sh` by your own setting.
By moving toolchain out or rename it ,one can make toolchain independent
from ABACUS repo, make dependencies package more independent and flexible.
2. Users can manually change `pkg_install_dir` variable 
in `scripts/stage*/install*` to change the installation directory 
of each packages, which may let the installation more fiexible.
3. Users can manually change `INSTALL` variable in `scripts/common_vars.sh`
to change the installation directory of all packages, which may let the
installation more fiexible.


## More
More infomation can be read from `Details.md`, 
which is merely easily refined from cp2k-toolchain README.