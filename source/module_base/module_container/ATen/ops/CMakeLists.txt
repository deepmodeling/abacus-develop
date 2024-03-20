file(GLOB ATen_OPS_SRCS "*.cpp")

set(ATen_CPU_SRCS ${ATen_CPU_SRCS} ${ATen_OPS_SRCS} PARENT_SCOPE)

if(BUILD_TESTING)
    add_subdirectory(test)
endif()