# =============================================================================
# Generate Build Information Header (Including Information Collection)
# ==============================================================================

# include_guard(GLOBAL)

function(setup_build_info)
    message(STATUS "Setting up build information...")

    include(cmake/CollectBuildInfoVars.cmake)

set(BUILD_INFO_TEMPLATE "${CMAKE_SOURCE_DIR}/source/source_io/build_info.h.in")
set(BUILD_INFO_OUTPUT   "${CMAKE_BINARY_DIR}/source/source_io/build_info.h")

configure_file(
    ${BUILD_INFO_TEMPLATE}
    ${BUILD_INFO_OUTPUT}
    @ONLY
)

    # add_library(BuildInfo::Headers INTERFACE IMPORTED GLOBAL)
    # target_include_directories(BuildInfo::Headers
    #     INTERFACE
    #         ${CMAKE_BINARY_DIR}/source/source_io
    # )

    message(STATUS "Build info header configured: ${BUILD_INFO_OUTPUT}")
endfunction()
