#!/bin/bash

# 脚本用于修复 get_device_type 函数调用，移除不必要的 ctx 参数

# 搜索所有调用 get_device_type 的文件
files=$(grep -r "get_device_type" --include="*.cpp" --include="*.hpp" --include="*.h" source/)

# 遍历每个文件
while IFS= read -r line; do
    # 提取文件路径
    file=$(echo "$line" | cut -d: -f1)
    
    # 打印正在处理的文件
    echo "Processing: $file"
    
    # 替换 get_device_type<Device>(ctx) 为 get_device_type<Device>()
    sed -i 's/get_device_type<\([^>]*\)>\(\s*\)(\s*\([^)]*\)\s*)/get_device_type<\1>\2()/g' "$file"
done <<< "$files"

echo "Fix completed!"
