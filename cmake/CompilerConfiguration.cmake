include(CheckCXXCompilerFlag)

# Match the existing fallback when no build type is specified, but keep the
# options scoped to ABACUS targets.
if(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
  target_compile_options(abacus_compile_requirements INTERFACE -O3 -g)
endif()

target_compile_options(abacus_compile_requirements INTERFACE
  "$<$<COMPILE_LANG_AND_ID:CXX,Intel,IntelLLVM>:-fp-model=strict;-Wno-write-strings>")

if(ENABLE_NATIVE_OPTIMIZATION AND NOT CMAKE_CROSSCOMPILING)
  set(_abacus_native_candidates)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
     OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(_abacus_native_candidates -march=native -mcpu=native)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(_abacus_native_candidates -xHOST)
  endif()
  foreach(_flag IN LISTS _abacus_native_candidates)
    string(MAKE_C_IDENTIFIER "ABACUS_HAS_${_flag}" _var)
    check_cxx_compiler_flag("${_flag}" ${_var})
    if(${_var})
      target_compile_options(abacus_compile_requirements INTERFACE
        "$<$<COMPILE_LANGUAGE:CXX>:${_flag}>")
      break()
    endif()
  endforeach()
endif()

if(USE_CUDA)
  target_compile_options(abacus_compile_requirements INTERFACE
    "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:-g;-G>")
  if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    separate_arguments(_abacus_openmp_cxx_flags NATIVE_COMMAND
                       "${OpenMP_CXX_FLAGS}")
    foreach(_flag IN LISTS _abacus_openmp_cxx_flags)
      target_compile_options(abacus_compile_requirements INTERFACE
        "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=${_flag}>")
    endforeach()
  endif()
endif()

unset(_abacus_native_candidates)
unset(_abacus_openmp_cxx_flags)
unset(_flag)
unset(_var)
