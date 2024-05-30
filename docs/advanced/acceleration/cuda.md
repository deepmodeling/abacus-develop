# CUDA GPU Implementations

In ABACUS, we provide the option to use GPU devices to accelerate performance. The implementation of GPU acceleration differs between PW basis and LCAO basis. Specifically, under PW basis, it has the following features:

- **Full gpu implementations**: During the SCF progress, `Psi`, `Hamilt`, `Hsolver`, `DiagCG`, and `DiagoDavid` classes are stored or calculated by the GPU devices.

- **Electronic state data**: (e.g. electronic density) are moved from the GPU to the CPU(s) every scf step.

- **Acclerated by the NVIDIA libraries**: `cuBLAS` for common linear algebra calculations, `cuSolver` for eigen values/vectors, and `cuFFT` for the conversions between the real and recip spaces.

- **Multi GPU supprted**: Using multiple MPI tasks will often give the best performance. Note each MPI task will be bind to a GPU device with automatically computing load balancing.

- **Parallel strategy**: K point parallel.

Unlike PW basis, only the grid integration module (module_gint) and the diagonalization of the Hamiltonian matrix (module_hsolver) have been implemented with GPU acceleration under LCAO basis, and the acceleration is limited to gamma only calculation. Additionally, LCAO basis does not support multi-GPU acceleration. Both the grid integration module and the Hamiltonian matrix solver only support acceleration on a single GPU.

## Required hardware/software

To compile and use ABACUS in CUDA mode, you currently need to have an NVIDIA GPU and install the corresponding NVIDIA CUDA toolkit software on your system (this is only tested on Linux and unsupported on Windows):

- Check if you have an NVIDIA GPU: cat /proc/driver/nvidia/gpus/*/information

- Go to https://developer.nvidia.com/cuda-downloads

- Install a driver and toolkit appropriate for your system (SDK is not necessary)


## Building ABACUS with the GPU support:

Check the [Advanced Installation Options](https://abacus-rtd.readthedocs.io/en/latest/advanced/install.html#build-with-cuda-support) for the installation of CUDA version support.

## Run with the GPU support by editing the INPUT script:

In `INPUT` file we need to set the value keyword [device](../input_files/input-main.md#device) to be `gpu`.

## Examples
We provides [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/gpu) of gpu calculations.

## Known limitations
PW basis:
- CG, BPCG and Davidson methods are supported, so the input keyword `ks_solver` can take the values `cg`, `bpcg` or `dav`.
- Only k point parallelization is supported, so the input keyword `kpar` will be set to match the number of MPI tasks automatically.
- By default, CUDA architectures 60, 70, 75, 80, 86, and 89 are compiled (if supported). It can be overriden using the CMake variable [`CMAKE_CUDA_ARCHITECTURES`](https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html) or the environmental variable [`CUDAARCHS`](https://cmake.org/cmake/help/latest/envvar/CUDAARCHS.html).

LCAO basis:
- Does not support multi-k calculation, so if the input keyword `device` is set to `gpu`, the input keyword `gamma_only` can only take the value `1`.
- Does not support multi-GPU acceleration.
