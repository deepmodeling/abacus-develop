if(ENABLE_FLOAT_FFTW)
  list(APPEND FFT_SRC ../../source_base/module_fft/fft_cpu_float.cpp)
endif()
if(USE_CUDA)
  list(APPEND FFT_SRC ../../source_base/module_fft/fft_cuda.cpp)
endif()
if(USE_ROCM)
  list(APPEND FFT_SRC ../../source_base/module_fft/fft_rocm.cpp)
endif()
if(USE_DSP)
  list(APPEND FFT_SRC ../../source_base/module_fft/fft_dsp.cpp
       ../../source_base/module_fft/fft_dsp_float.cpp pw_transform_k_dsp.cpp)
endif()

list(
  APPEND
  objects
  pw_basis.cpp
  pw_grid_geometry.cpp
  pw_basis_k.cpp
  pw_basis_sup.cpp
  pw_distributeg.cpp
  pw_distg_method1.cpp
  pw_distg_method2.cpp
  pw_distributer.cpp
  pw_init.cpp
  pw_transform.cpp
  pw_transform_gpu.cpp
  pw_transform_k.cpp
  ../../source_base/module_fft/fft_bundle.cpp
  ../../source_base/module_fft/fft_cpu.cpp
  ${FFT_SRC})

add_library(planewave OBJECT ${objects})

if(USE_DSP)
  target_compile_definitions(
    planewave PUBLIC FFT_DAT_DIR="${MTBLAS_FFT_DIR}/datfile/mt_fft_blas.dat")
endif()
if(USE_CUFFTDX)
  set(fused_sources kernels/cuda/fused_veff.cu)
  set(fused_declarations "")
  set(fused_cases "")
  set(FFT_MIN_LENGTH 16)
  set(FFT_MAX_LENGTH 256)
  set(MIN_FUSED_FFT_SM 70) # SM 7.0
  set(SM_VERSION_ENCODING_FACTOR 10)
  set(fused_architectures ${CMAKE_CUDA_ARCHITECTURES})
  if("native" IN_LIST fused_architectures)
    set(fused_architectures ${CMAKE_CUDA_ARCHITECTURES_NATIVE})
  elseif("all" IN_LIST fused_architectures)
    set(fused_architectures ${CMAKE_CUDA_ARCHITECTURES_ALL})
  elseif("all-major" IN_LIST fused_architectures)
    set(fused_architectures ${CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR})
  endif()
  set(fused_sms "")
  foreach(architecture IN LISTS fused_architectures)
    if(NOT architecture MATCHES "^([0-9]+)([af]?)(-(real|virtual))?$")
      message(FATAL_ERROR "USE_CUFFTDX requires explicit, native, all, or all-major CUDA architectures.")
    endif()
    set(sm ${CMAKE_MATCH_1})
    if(sm GREATER_EQUAL ${MIN_FUSED_FFT_SM})
      list(APPEND fused_sms ${sm})
      list(APPEND fused_targets_${sm} ${architecture})
    endif()
  endforeach()
  list(REMOVE_DUPLICATES fused_sms)
  list(SORT fused_sms COMPARE NATURAL ORDER DESCENDING)
  foreach(FFT_LENGTH RANGE ${FFT_MIN_LENGTH} ${FFT_MAX_LENGTH})
    set(remainder ${FFT_LENGTH})
    foreach(factor 2 3 5)
      math(EXPR residue "${remainder} % ${factor}")
      while(residue EQUAL 0)
        math(EXPR remainder "${remainder} / ${factor}")
        math(EXPR residue "${remainder} % ${factor}")
      endwhile()
    endforeach()
    if(remainder EQUAL 1)
      string(APPEND fused_cases "        case ${FFT_LENGTH}:\n")
      foreach(FFT_ARCH IN LISTS fused_sms)
        math(EXPR FFT_SM "${FFT_ARCH} * ${SM_VERSION_ENCODING_FACTOR}")
        set(generated "${CMAKE_CURRENT_BINARY_DIR}/fused_fft_${FFT_LENGTH}_${FFT_ARCH}.cu")
        configure_file(kernels/cuda/fused_fft_size.cu.in ${generated} @ONLY)
        set(target planewave_fft_${FFT_LENGTH}_${FFT_ARCH})
        add_library(${target} OBJECT ${generated})
        target_link_libraries(${target} PRIVATE mathdx::cufftdx)
        set_target_properties(${target} PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
                              CUDA_ARCHITECTURES "${fused_targets_${FFT_ARCH}}")
        target_sources(${ABACUS_BIN_NAME} PRIVATE $<TARGET_OBJECTS:${target}>)
        string(APPEND fused_declarations
          "const FusedFFT* fused_fft_${FFT_LENGTH}_${FFT_ARCH}(const cudaDeviceProp&);\n")
        string(APPEND fused_cases
          "            if (sm >= ${FFT_SM}) {\n"
          "                if (const auto* fft = fused_fft_${FFT_LENGTH}_${FFT_ARCH}(properties)) { return fft; }\n"
          "            }\n")
      endforeach()
      string(APPEND fused_cases "            return nullptr;\n")
    endif()
  endforeach()
  set(registry "${CMAKE_CURRENT_BINARY_DIR}/fused_fft_registry.cpp")
  configure_file(kernels/cuda/fused_fft_registry.cpp.in ${registry} @ONLY)
  list(APPEND fused_sources ${registry})
  add_library(planewave_fused_fft OBJECT ${fused_sources})
  target_link_libraries(planewave_fused_fft PRIVATE mathdx::cufftdx)
  set_target_properties(planewave_fused_fft PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
  target_sources(${ABACUS_BIN_NAME} PRIVATE $<TARGET_OBJECTS:planewave_fused_fft>)
endif()

if(ENABLE_COVERAGE)
  add_coverage(planewave)
endif()

if(BUILD_TESTING)
  add_subdirectory(test_serial)
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(kernels/test)
    add_subdirectory(test_gpu)
  endif()
endif()
