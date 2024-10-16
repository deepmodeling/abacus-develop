#ifndef DIAG_CONST_NUMS
#define DIAG_CONST_NUMS
#include "module_base/module_device/memory_op.h"

template <typename T>
struct const_nums
{
    const_nums();
    base_device::DEVICE_CPU* cpu_ctx = {};
    T* zero;
    T* one;
    T* neg_one;
};

#endif