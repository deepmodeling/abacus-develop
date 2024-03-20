set(BASE_CORE_CPU_SRCS
    base/core/refcount.cpp
    base/core/cpu_allocator.cpp
)

if (USE_CUDA OR USE_ROCM)
    list(APPEND BASE_CORE_CPU_SRCS 
        base/core/gpu_allocator.cpp
    )
endif()

set(ATen_CPU_SRCS ${ATen_CPU_SRCS} ${BASE_CORE_CPU_SRCS} PARENT_SCOPE)