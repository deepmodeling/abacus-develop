# CUDA GPU Implementations

In ABACUS, we provide the option to use GPU devices to accelerate performance. The implementation of GPU acceleration differs between PW basis and LCAO basis. Specifically, under PW basis, it has the following features:

- **Full gpu implementations**: During the SCF progress, `Psi`, `Hamilt`, `Hsolver`, `DiagCG`, and `DiagoDavid` classes are stored or calculated by the GPU devices.

- **Electronic state data**: (e.g. electronic density) are moved from the GPU to the CPU(s) every scf step.

- **Accelerated by the NVIDIA libraries**: `cuBLAS` for common linear algebra calculations, `cuSolver` for eigen values/vectors, and `cuFFT` for the conversions between the real and recip spaces.

- **Multi GPU supprted**: Using multiple MPI tasks will often give the best performance. Note each MPI task will be bind to a GPU device with automatically computing load balancing.

- **Parallel strategy**: K point parallel.

Unlike PW basis, only the grid integration module (module_gint) and the diagonalization of the Hamiltonian matrix (module_hsolver) have been implemented with GPU acceleration under LCAO basis.

## Required hardware/software

To compile and use ABACUS in CUDA mode, you currently need to have an NVIDIA GPU and install the corresponding NVIDIA CUDA toolkit software on your system (this is only tested on Linux and unsupported on Windows):

- Check if you have an NVIDIA GPU: cat /proc/driver/nvidia/gpus/*/information

- Go to https://developer.nvidia.com/cuda-downloads

- Install a driver and toolkit appropriate for your system (SDK is not necessary)


## Building ABACUS with the GPU support:

Check the [Advanced Installation Options](https://abacus-rtd.readthedocs.io/en/latest/advanced/install.html#build-with-cuda-support) for the installation of CUDA version support.

Setting both USE_ELPA and USE_CUDA to ON does not automatically enable ELPA to run on GPUs. ELPA support for GPUs needs to be enabled when ELPA is compiled. [enable GPU support](https://github.com/marekandreas/elpa/blob/master/documentation/INSTALL.md).

The ABACUS program will automatically determine whether the current ELPA supports GPU based on the elpa/elpa_configured_options.h header file. Users can also check this header file to determine the GPU support of ELPA in their environment. ELPA introduced a new API elpa_setup_gpu in version 2023.11.001. So if you want to enable ELPA GPU in ABACUS, the ELPA version must be greater than or equal to 2023.11.001.

## Run with the GPU support by editing the INPUT script:

In `INPUT` file we need to set the input parameter [device](../input_files/input-main.md#device) to `gpu`. If this parameter is not set, ABACUS will try to determine if there are available GPUs.
- Set `ks_solver`: For the PW basis, CG, BPCG and Davidson methods are supported on GPU; set the input parameter [ks_solver](../input_files/input-main.md#ks_solver) to `cg`, `bpcg` or `dav`. For the LCAO basis, `cusolver` and `elpa` is supported on GPU.
- **multi-card**: ABACUS allows for multi-GPU acceleration. If you have multiple GPU cards, you can run ABACUS with several MPI processes, and each process will utilize one GPU card. For example, the command `mpirun -n 2 abacus` will by default launch two GPUs for computation. If you only have one card, this command will only start one GPU. 

## Examples
We provides [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/gpu) of gpu calculations.

## Known limitations
PW basis:
- Only k point parallelization is supported, so the input keyword `kpar` will be set to match the number of MPI tasks automatically.
- By default, CUDA architectures 60, 70, 75, 80, 86, and 89 are compiled (if supported). It can be overriden using the CMake variable [`CMAKE_CUDA_ARCHITECTURES`](https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html) or the environmental variable [`CUDAARCHS`](https://cmake.org/cmake/help/latest/envvar/CUDAARCHS.html).
LCAO basis:
- Unless there is a specific reason, avoid using multiple GPUs, as it can be slower than using a single GPU. This is because the generalized eigenvalue solution of the LCAO basis set will incur additional communication overhead when calculated on multiple cards. When the memory limit of a GPU card makes it insufficient to complete the task, it is recommended to use multiple cards for calculation.
- When using elpa on GPUs, some ELPA internal logs will be output.