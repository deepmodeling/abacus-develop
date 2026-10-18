###############################################################################
# - Find LibRI
# Find the native LibRI files.
#
#  LIBRI_FOUND - True if LibRI is found.
#  LIBRI_DIR - Where to find LibRI files.
#  LibRI::LibRI - Header-only target, including optional DDLA requirements.

find_path(LIBRI_DIR
    include/RI/version.h
    HINTS ${LIBRI_DIR}
    HINTS ${LibRI_DIR}
    HINTS ${libri_DIR}
)

if(NOT LIBRI_DIR)
    if(ENABLE_LIBRI_GPU)
        message(FATAL_ERROR
            "ENABLE_LIBRI_GPU requires an explicitly supplied, patched LIBRI_DIR; "
            "see docs/advanced/install.md")
    endif()
    include(FetchContent)
    FetchContent_Declare(
        LibRI
        URL https://codeload.github.com/abacusmodeling/LibRI/tar.gz/b0eff7a
    )
    FetchContent_Populate(LibRI)
    set(LIBRI_DIR ${libri_SOURCE_DIR})
endif()
# Handle the QUIET and REQUIRED arguments and
# set LIBRI_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibRI DEFAULT_MSG LIBRI_DIR)

if(LibRI_FOUND AND NOT TARGET LibRI::LibRI)
    find_package(LibComm REQUIRED)
    find_package(cereal REQUIRED CONFIG)

    add_library(LibRI::LibRI INTERFACE IMPORTED)
    set_target_properties(LibRI::LibRI PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBRI_DIR}/include;${LIBCOMM_DIR}/include"
        INTERFACE_COMPILE_FEATURES cxx_std_14
        INTERFACE_LINK_LIBRARIES "MPI::MPI_CXX;cereal::cereal")
    if(ENABLE_OPENMP)
        target_link_libraries(LibRI::LibRI INTERFACE OpenMP::OpenMP_CXX)
    endif()
    # Reuse ABACUS's selected math backend. Special platforms without these
    # targets retain their existing numerical libraries in the final link.
    foreach(_libri_math_target IN ITEMS LAPACK::LAPACK BLAS::BLAS)
        if(TARGET ${_libri_math_target})
            target_link_libraries(LibRI::LibRI INTERFACE ${_libri_math_target})
        endif()
    endforeach()

    if(ENABLE_LIBRI_GPU)
        find_package(DDLA REQUIRED)
        target_compile_definitions(LibRI::LibRI INTERFACE __GPU_RI __DDLA_RI)
        target_compile_features(LibRI::LibRI INTERFACE cxx_std_17)
        target_link_libraries(LibRI::LibRI INTERFACE DDLA::DDLA)

        include(CMakePushCheckState)
        include(CheckCXXSourceCompiles)
        cmake_push_check_state(RESET)
        set(CMAKE_REQUIRED_LIBRARIES LibRI::LibRI)
        # Recheck when users change dependency paths or apply compatibility fixes.
        unset(LIBRI_GPU_API_WORKS CACHE)
        check_cxx_source_compiles([=[
            #include <complex>
            #include <RI/ri/LRI.h>
            template <typename T> void check_gpu()
            {
                RI::LRI<int, int, 3, T> lri;
                lri.cal_mode = RI::LRI_Cal_Mode::GPU;
                std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<T>>> result;
                lri.cal_loop3({}, result, 1.0);
            }
            int main()
            {
                check_gpu<double>();
                check_gpu<std::complex<double>>();
            }
        ]=] LIBRI_GPU_API_WORKS)
        cmake_pop_check_state()
        if(NOT LIBRI_GPU_API_WORKS)
            message(FATAL_ERROR
                "LibRI/DDLA GPU compile/link check failed. Supply compatible, patched "
                "LibRI and CUDA DDLA installations; see docs/advanced/install.md "
                "and the CMake configure log.")
        endif()
    endif()
endif()

mark_as_advanced(LIBRI_DIR)
