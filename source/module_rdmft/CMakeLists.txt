if(ENABLE_LCAO)
    add_library(
        rdmft
        OBJECT
        rdmft.cpp
        rdmft_tools.cpp
    )
endif()

# if(ENABLE_COVERAGE)
#   add_coverage(psi)
#   add_coverage(psi_initializer)
# endif()

# if (BUILD_TESTING)
#   add_subdirectory(kernels/test)
#   add_subdirectory(test)
# endif()
