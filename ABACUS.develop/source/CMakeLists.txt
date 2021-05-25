add_subdirectory(src_global)
add_subdirectory(module_ORB)

add_executable(main module_ORB/main.cpp)
target_link_options(main PRIVATE -mkl)
target_link_libraries(main PRIVATE orb global)
