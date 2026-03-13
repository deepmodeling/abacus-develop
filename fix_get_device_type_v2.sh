#!/bin/bash

# 脚本用于修复 get_device_type 函数调用，去掉模板参数

# 定义要修改的文件列表
files=(
    "source/source_pw/module_stodft/sto_forces.cpp"
    "source/source_pw/module_stodft/sto_wf.cpp"
    "source/source_pw/module_pwdft/fs_nonlocal_tools.cpp"
    "source/source_pw/module_pwdft/onsite_projector.cpp"
    "source/source_pw/module_pwdft/onsite_proj_tools.cpp"
    "source/source_pw/module_pwdft/nonlocal_maths.hpp"
    "source/source_pw/module_pwdft/op_pw_ekin.cpp"
    "source/source_pw/module_pwdft/stress_loc.cpp"
    "source/source_pw/module_pwdft/stress_cc.cpp"
    "source/source_hsolver/test/hsolver_pw_sup.h"
    "source/source_pw/module_pwdft/fs_kin_tools.cpp"
    "source/source_pw/module_pwdft/forces_scc.cpp"
    "source/source_pw/module_pwdft/forces_cc.cpp"
    "source/source_pw/module_pwdft/forces.cpp"
    "source/source_hsolver/diago_iter_assist.cpp"
    "source/source_hsolver/diago_david.cpp"
    "source/source_hsolver/diago_dav_subspace.cpp"
    "source/source_esolver/esolver_ks_pw.cpp"
    "source/source_base/math_chebyshev.cpp"
    "source/source_pw/module_pwdft/structure_factor_k.cpp"
    "source/source_base/module_device/test/device_test.cpp"
)

# 遍历每个文件
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "Processing: $file"
        # 替换 get_device_type<Device>( 为 get_device_type(
        sed -i 's/get_device_type<Device>(/get_device_type(/g' "$file"
        # 替换 get_device_type<base_device::DEVICE_CPU>( 为 get_device_type(
        sed -i 's/get_device_type<base_device::DEVICE_CPU>(/get_device_type(/g' "$file"
        # 替换 get_device_type<base_device::DEVICE_GPU>( 为 get_device_type(
        sed -i 's/get_device_type<base_device::DEVICE_GPU>(/get_device_type(/g' "$file"
    else
        echo "File not found: $file"
    fi
done

echo "Fix completed!"
