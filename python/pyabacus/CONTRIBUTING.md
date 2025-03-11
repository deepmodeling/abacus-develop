# Developer Guide

## Introduction

Welcome to the `pyabacus` project! This document provides guidelines and instructions for developers who want to contribute to this project.

`pyabacus` is a Python interface for the ABACUS package. It provides a high-level Python API for interacting with the ABACUS library, allowing users to perform electronic structure calculations and analyze the results using Python.

<!-- toc -->

- [Developer Guide](#developer-guide)
  - [Introduction](#introduction)
  - [Project Structure](#project-structure)
    - [Root CMake Configuration](#root-cmake-configuration)
    - [Module CMake Configuration](#module-cmake-configuration)
  - [Development Process](#development-process)
  - [Conclusion](#conclusion)

<!-- tocstop -->

**If you are new to the project**, please refer to the [README.md](./README.md) file for an overview of the project and its goals.

**If you are already familiar with the project and want to contribute**, this guide will help you understand the project structure, development process, and best practices for contributing code.

**If you have any questions or need help**, feel free to reach out to the maintainers or create an issue in the repository.

**Please feel free to contribute to this guide** by submitting a pull request with any improvements or additional information.

Let's get started!

## Project Structure

The project is organized as follows:

```
pyabacus/
├── CMakeLists.txt
└── src
    ├── pyabacus
    │   └── {your_module}
    │       ├── {interface}.py
    │       └── __init__.py
    └── {your_module}
        ├── {your_code}.cpp
        └── CMakeLists.txt
```

Our project is built using [pybind11](http://github.com/pybind/pybind11) and [scikit-build-core](https://scikit-build-core.readthedocs.io/) for facilitating the `CMake` build toolchain. So the `CMakeLists.txt` configuration is the key to thoroughly understanding the project structure.

### Root CMake Configuration

The `CMakeLists.txt` in root directory is the main configuration file for the pyabacus project. It sets up the project, finds necessary dependencies, configures build options, and includes subdirectories for different modules. Below is a detailed explanation of each section of the file:

```cmake
cmake_minimum_required(VERSION 3.15...3.26)

# Project settings 
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)
```
- This section sets the project name, version, and the programming languages used (C++ in this case). The project name and version are obtained from the `SKBUILD_PROJECT_NAME` and `SKBUILD_PROJECT_VERSION` variables, respectively.

```cmake
# Find Python and pybind11
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
```
- This section finds the required Python and pybind11 packages. 

```cmake
# Set source path
set(ABACUS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../source")
set(BASE_PATH "${ABACUS_SOURCE_DIR}/module_base")
set(NAO_PATH "${ABACUS_SOURCE_DIR}/module_basis/module_nao")
set(HSOLVER_PATH "${ABACUS_SOURCE_DIR}/module_hsolver")
set(PSI_PATH "${ABACUS_SOURCE_DIR}/module_psi")
set(ENABLE_LCAO ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")
```
- This section sets various source paths and configuration options. It defines the paths to different modules and appends the custom CMake module path.

```cmake
# Add math_libs 
if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
endif()
if(MKLROOT)
  set(MKL_INTERFACE lp64)
  set(ENABLE_MPI ON)
  if (ENABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
  endif()

  set(USE_OPENMP ON)
  if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_link_options(${OpenMP_CXX_LIBRARIES})
  endif()
  find_package(MKL REQUIRED)
  add_definitions(-D__MKL)
  include_directories(${MKL_INCLUDE} ${MKL_INCLUDE}/fftw)

  if(NOT ENABLE_DEEPKS)
    list(APPEND math_libs IntelMKL::MKL)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    list(APPEND math_libs -lifcore)
  endif()
else()
    find_package(FFTW3 REQUIRED)
    add_compile_definitions(__FFTW3)
    find_package(LAPACK REQUIRED)
    include_directories(${FFTW3_INCLUDE_DIRS})
    list(APPEND math_libs FFTW3::FFTW3 LAPACK::LAPACK)

  if(ENABLE_LCAO)
    find_package(ScaLAPACK REQUIRED)
    list(APPEND math_libs ScaLAPACK::ScaLAPACK)
  endif()
endif()
```
- This section configures the math libraries. It checks for the presence of the Intel Math Kernel Library (MKL) and configures it if available. If MKL is not available, it falls back to using FFTW3 and LAPACK. It also configures MPI and OpenMP if enabled.

```cmake
# Add include directories
include_directories(
    ${BASE_PATH} 
    ${ABACUS_SOURCE_DIR}
    ${ABACUS_SOURCE_DIR}/module_base/module_container
    )
```
- This section adds the necessary include directories for the project.

```cmake
# Add basic libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Add base
set(BASE_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/base")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_base ${BASE_BINARY_DIR})
# Add parameter
set(PARAMETER_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/parameter")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_parameter ${PARAMETER_BINARY_DIR})
# Add orb
set(ORB_BINARY_DIR "${PROJECT_SOURCE_DIR}/build/orb")
add_subdirectory(${ABACUS_SOURCE_DIR}/module_basis/module_ao ${ORB_BINARY_DIR})
```
- This section sets the position-independent code flag and adds subdirectories for the base, parameter, and orb modules. It specifies the build directories for these modules.

```cmake
# Set RPATH
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
```
- This section sets the runtime search path (RPATH) for the Python site-packages directory. It uses a Python command to get the site-packages path and stores it in the `PYTHON_SITE_PACKAGES` variable.

```cmake
# Set package name to pyabacus
set(TARGET_PACK pyabacus)
set(CMAKE_INSTALL_RPATH "${PYTHON_SITE_PACKAGES}/${TARGET_PACK}")
```
- This section sets the package name to `pyabacus` and configures the install RPATH to include the Python site-packages directory.

```cmake
# Add subdirectories for submodules
add_subdirectory(${PROJECT_SOURCE_DIR}/src/hsolver)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ModuleBase)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ModuleNAO)
```
- This section adds subdirectories for modules. Each subdirectory contains its own `CMakeLists.txt` for further configuration.

By following this structure, the `CMakeLists.txt` file ensures that all necessary dependencies are found, configured, and included in the build process. It also sets up the project environment and includes submodules for different components of the `pyabacus` project.

### Module CMake Configuration

I'll show you a `CMakeLists.txt` for example (pyabacus.hsolver)

```cmake
# Add diago shared library
list(APPEND _diago
    ${HSOLVER_PATH}/diago_dav_subspace.cpp
    ${HSOLVER_PATH}/diago_david.cpp
    ${HSOLVER_PATH}/diag_const_nums.cpp
    ${HSOLVER_PATH}/diago_iter_assist.cpp
    ${HSOLVER_PATH}/kernels/dngvd_op.cpp
    ${BASE_PATH}/kernels/math_kernel_op.cpp
    ${BASE_PATH}/kernels/math_op.cpp
    ${BASE_PATH}/module_device/device.cpp
    ${BASE_PATH}/module_device/memory_op.cpp
    ${PSI_PATH}/psi.cpp
)
add_library(diagopack SHARED ${_diago})
target_link_libraries(diagopack
    base
    parameter
    container
    orb
    ${math_libs}
    ${OpenBLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
)

list(APPEND pymodule_hsolver
    ${PROJECT_SOURCE_DIR}/src/hsolver/py_hsolver.cpp
)

# Use pybind11 to add python module
pybind11_add_module(_hsolver_pack MODULE ${pymodule_hsolver})
# Link your dependencies and pybind11 libraries to your module 
target_link_libraries(_hsolver_pack PRIVATE pybind11::headers diagopack)
target_compile_definitions(_hsolver_pack PRIVATE VERSION_INFO=${PROJECT_VERSION})

set_target_properties(diagopack PROPERTIES INSTALL_RPATH "$ORIGIN")
set_target_properties(_hsolver_pack PROPERTIES INSTALL_RPATH "$ORIGIN")

# Install your module package to destination path
install(TARGETS _hsolver_pack diagopack DESTINATION ${TARGET_PACK}/hsolver)
```

You can refer to the `CMakeLists.txt` files in other modules for guidance on how to configure your module.

## Development Process

To contribute to the `pyabacus` project, follow these steps:

1.  **Check the issues**:
    - Look for issues to ensure that you are not working on something that is already in progress.
    - If you want to work on a new feature or bug fix, create an issue first to discuss it with the maintainers.

2. **Create a new folder for your module**:
   - If you want to add a new module with pure Python code, create a new folder in the `src/pyabacus` directory.
    - If you want to add a new module with C++ code, create a new folder in the `src` directory and a corresponding directory in the `src/pyabacus` directory.

3. **Write source code using pybind11**:
   - Follow the structure of other modules.
   - Manage dependencies and installation paths in the `CMakeLists.txt` file.

3. **Modify `src/pyabacus/__init__.py`**:
   - Add the name of your module to the `__submodules__` list and import the module in the `__getattr__` function.

   ```python
   from __future__ import annotations

   __submodules__ = ["ModuleBase", "ModuleNAO", "hsolver", "{module_name}"]

   __all__ = list(__submodules__)

   def __getattr__(attr):
       if attr == "ModuleBase":
           import pyabacus.ModuleBase as ModuleBase
           return ModuleBase
       elif attr == "ModuleNAO":
           import pyabacus.ModuleNAO as ModuleNAO
           return ModuleNAO
       elif attr == "hsolver":
           import pyabacus.hsolver as hsolver
           return hsolver
       elif attr == '{module_name}':
           import pyabacus.{module_name} as {module_name}
           return {module_name}
       else:
           raise AttributeError(f"module {__name__} has no attribute {attr}")
   ```

4. **Create two files in `src/pyabacus/{module_name}`**:
   - `__init__.py`: This file allows Python to recognize the folder as a module.
   - `_{module_name}.py`: This file is responsible for designing the Python interface (frontend).

   **Example `__init__.py`**:

   ```python
   from __future__ import annotations
   from ._{module_name} import *

   __all__ = ["{class_name}", "{func_name}", ...]
   ```

   **Example `_{module_name}.py`**:

   ```python
   from .{module_library_name} import {your_class} as _your_class, ...

   """
   Your class should inherit from the corresponding class in the C++ library.
   All methods should be overridden to provide type hints and auto-completion.
   You can use the `super()` method to call the base class(C++ class) methods.
   """
   class {your_class}(_your_class):
       def __init__(self) -> None:
           super().__init__()
       
       def foo(self, arg1, arg2, ...) -> RetType:
           return super().foo(arg1, arg2, ...)
       
       def bar(self, arg1, arg2, ...):
           super().bar(arg1, arg2, ...)
   ```

   For a class, if you do not declare the interface in the frontend, the IDE will not provide type hints and auto-completion. However, if the interface name matches the name binding in pybind11, it will be overridden. To address this, you can use the method as shown above.

5. **Handle overloaded functions in C++**:
   - Since Python does not support function overloading with different parameters, use the following method:

   ```python
   @overload
   def foo(self, x: float) -> float: ...
   @overload
   def foo(self, n: int, x: float, y: float) -> float: ...

   def foo(self, *args, **kwargs):
       return super().foo(*args, **kwargs)
   ```

**Example Python Interface**:

   ```python
   class diag_comm_info(_diag_comm_info):
       def __init__(self, rank: int, nproc: int):
           super().__init__(rank, nproc)
       
       @property
       def rank(self) -> int:
           return super().rank
       
       @property
       def nproc(self) -> int:
           return super().nproc

   class Sphbes(_Sphbes):
       def __init__(self) -> None: 
           super().__init__()
           
       @overload
       @staticmethod
       def sphbesj(l: int, x: float) -> float: ...
       @overload
       @staticmethod
       def sphbesj(
           n: int, 
           r: NDArray[np.float64], 
           q: int, 
           l: int, 
           jl: NDArray[np.float64]
       ) -> None: ...
       
       def sphbesj(self, *args, **kwargs): 
           return super().sphbesj(*args, **kwargs)
           
       @overload
       @staticmethod
       def dsphbesj(l: int, x: float) -> float: ...
       @overload
       @staticmethod
       def dsphbesj(
           n: int, 
           r: NDArray[np.float64], 
           q: int, 
           l: int, 
           djl: NDArray[np.float64]
       ) -> None: ...
       
       def dsphbesj(self, *args, **kwargs):
           return super().dsphbesj(*args, **kwargs)
           
       @staticmethod
       def sphbes_zeros(l: int, n: int, zeros: NDArray[np.float64]) -> None: 
           super().sphbes_zeros(l, n, zeros)
   ```

## Conclusion

By following this guide, you can effectively contribute to the `pyabacus` project. Ensure that you follow the structure and conventions outlined here to maintain consistency and readability in the codebase. Happy coding!
