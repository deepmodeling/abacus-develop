if (ENABLE_FLOAT_FFTW)
  list (APPEND FFT_SRC
    module_fft/fft_cpu_float.cpp
  )
endif()
if (USE_CUDA)
  list (APPEND FFT_SRC
    module_fft/fft_cuda.cpp
  )
endif()
if (USE_ROCM)
  list (APPEND FFT_SRC
    module_fft/fft_rocm.cpp
  )
endif()

list(APPEND objects
    pw_basis.cpp
    pw_basis_k.cpp
    pw_basis_sup.cpp
    pw_distributeg.cpp
    pw_distributeg_method1.cpp
    pw_distributeg_method2.cpp
    pw_distributer.cpp
    pw_init.cpp
    pw_transform.cpp
    pw_transform_k.cpp
    module_fft/fft_bundle.cpp
    module_fft/fft_cpu.cpp
    ${FFT_SRC}
)

add_library(
    planewave
    OBJECT
    ${objects}
)

if(USE_CUFFTDX)
  set(fused_sources kernels/cuda/fused_veff.cu)
  set(fused_declarations "")
  set(fused_cases "")
  # Enumerate the smooth lengths selected by PW_Basis::initgrids. Each
  # translation unit instantiates the same kernels and can compile in parallel.
  foreach(FFT_LENGTH RANGE 16 256)
    set(remainder ${FFT_LENGTH})
    foreach(factor 2 3 5)
      math(EXPR residue "${remainder} % ${factor}")
      while(residue EQUAL 0)
        math(EXPR remainder "${remainder} / ${factor}")
        math(EXPR residue "${remainder} % ${factor}")
      endwhile()
    endforeach()
    if(remainder EQUAL 1)
      set(generated "${CMAKE_CURRENT_BINARY_DIR}/fused_fft_${FFT_LENGTH}.cu")
      configure_file(kernels/cuda/fused_fft_size.cu.in ${generated} @ONLY)
      list(APPEND fused_sources ${generated})
      string(APPEND fused_declarations "const FusedFFT* fused_fft_${FFT_LENGTH}();\n")
      string(APPEND fused_cases "        case ${FFT_LENGTH}: return fused_fft_${FFT_LENGTH}();\n")
    endif()
  endforeach()
  set(registry "${CMAKE_CURRENT_BINARY_DIR}/fused_fft_registry.cpp")
  file(WRITE ${registry}
    "#include \"module_basis/module_pw/kernels/fused_veff.h\"\n"
    "namespace ModulePW {\n${fused_declarations}"
    "const FusedFFT* fused_fft(int length) {\n    switch (length) {\n${fused_cases}"
    "        default: return nullptr;\n    }\n}\n}\n")
  list(APPEND fused_sources ${registry})
  add_library(planewave_fused_fft OBJECT ${fused_sources})
  target_link_libraries(planewave_fused_fft PRIVATE mathdx::cufftdx)
  set_target_properties(planewave_fused_fft PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES 70)
  target_sources(${ABACUS_BIN_NAME} PRIVATE $<TARGET_OBJECTS:planewave_fused_fft>)
endif()

if(ENABLE_COVERAGE)
  add_coverage(planewave)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(test_serial)
    add_subdirectory(kernels/test)
  endif()
endif()
