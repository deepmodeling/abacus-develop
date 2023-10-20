# CUDA GPU Implementations

In ABACUS, we provide the option to use the GPU devices to accelerate the performance.
And it has the following general features:

- **Full gpu implementations**: During the SCF progress, `Psi`, `Hamilt`, `Hsolver`, `DiagCG`, and `DiagoDavid` classes are stored or calculated by the GPU devices.

- **Electronic state data**: (e.g. electronic density) are moved from the GPU to the CPU(s) every scf step.

- **Acclerated by the NVIDIA libraries**: `cuBLAS` for common linear algebra calculations, `cuSolver` for eigen values/vectors, and `cuFFT` for the conversions between the real and recip spaces.

- **Multi GPU supprted**: Using multiple MPI tasks will often give the best performance. Note each MPI task will be bind to a GPU device with automatically computing load balancing.

- **Parallel strategy**: K point parallel.

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

- CG, BPCG and Davidson methods are supported, so the input keyword `ks_solver` can take the values `cg`, `bpcg` or `dav`,
- Only PW basis is supported, so the input keyword `basis_type` can only take the value `pw`,
- Only k point parallelization is supported, so the input keyword `kpar` will be set to match the number of MPI tasks automatically.
- Supported CUDA architectures:
  - 60 # P100, 1080ti
  - 70 # V100
  - 75 # T4
  - 80 # A100, 3090

## FAQ
```
Q: Does the GPU implementations support atomic orbital basis sets?
A: Currently no.
```
