remove_definitions(-D__MPI)

# Unit tests for module_math

if(BUILD_TESTING)
    # Test for spherical Bessel functions
    AddTest(
        TARGET test_math_sphbes
        LIBS module_math
        SOURCES math_sphbes_test.cpp
    )

    # Test for Chebyshev polynomials
    AddTest(
        TARGET test_math_chebyshev
        LIBS module_math
        SOURCES math_chebyshev_test.cpp
    )

    # Test for B-spline
    AddTest(
        TARGET test_math_bspline
        LIBS module_math
        SOURCES math_bspline_test.cpp
    )

    # Test for polynomial interpolation
    AddTest(
        TARGET test_math_polyint
        LIBS module_math
        SOURCES math_polyint_test.cpp
    )

    # Test for integral
    AddTest(
        TARGET test_math_integral
        LIBS module_math
        SOURCES math_integral_test.cpp
    )

    # Test for Ylm real
    AddTest(
        TARGET test_math_ylmreal
        LIBS module_math
        SOURCES math_ylmreal_test.cpp
    )

    # Test for mymath
    AddTest(
        TARGET test_mymath
        LIBS module_math
        SOURCES mymath_test.cpp
    )

    # Test for spherical Bessel functions performance
    AddTest(
        TARGET test_perf_sphbes
        LIBS module_math
        SOURCES perf_sphbes_test.cpp
    )

    # Test for Lebedev-Laikov grid
    AddTest(
        TARGET test_lebedev_laikov
        LIBS module_math
        SOURCES test_lebedev_laikov.cpp
    )
endif()
