#!/bin/bash

# 检查文件是否存在
if [ $# -ne 1 ] || [ ! -f "$1" ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

file="$1"
found=0

# 逐行读取文件内容
while IFS= read -r line || [[ -n "$line" ]]; do

    # 检查是否包含关键字“key”
    if [[ "$line" == *"hello from set_target_mag for nspin_ = 4"* ]]; then
        found=2
        continue
    fi

    # 如果已经找到“key”，捕捉后续两行
    if [ $found -gt 0 ]; then
        # 存储两行数据
        #vector lines
        vectors[$found]="$line"
        found=$((found - 1))

        # 当捕捉到两行后进行处理
        if [ $found -eq 0 ]; then
            # 提取矢量数据
            vector1="${vectors[2]}"
            vector2="${vectors[1]}"

            # 清理矢量数据格式（假设格式为 [x,y,z] 或 类似形式）
            vector1=$(echo "$vector1" | sed 's/^0-9.,-//g' | tr ',' ' ')
            vector2=$(echo "$vector2" | sed 's/^0-9.,-//g' | tr ',' ' ')

            # 检查矢量格式是否正确
            if [ $(echo "$vector1" | wc -w) -ne 3 ] || [ $(echo "$vector2" | wc -w) -ne 3 ]; then
                echo "Invalid vector format, skipping..."
                continue
            fi

            # 将矢量转换为数组
            IFS=" " read -r -a vec1 <<< "$vector1"
            IFS=" " read -r -a vec2 <<< "$vector2"
            #echo "vec1: ${vec1[@]}" 
            #echo "vec2: ${vec2[@]}"

            # 计算点积
            dot_product=$(echo "scale=10; (${vec1[0]} * ${vec2[0]}) + (${vec1[1]} * ${vec2[1]}) + (${vec1[2]} * ${vec2[2]})" | bc)

            # 计算矢量的模长
            mod_vec1=$(echo "scale=10; sqrt((${vec1[0]}^2) + (${vec1[1]}^2) + (${vec1[2]}^2))" | bc -l)
            mod_vec2=$(echo "scale=10; sqrt((${vec2[0]}^2) + (${vec2[1]}^2) + (${vec2[2]}^2))" | bc -l)

            # 避免除以零的情况
            #if [ $(echo "$mod_vec1 * $mod_vec2" | bc) -le 0.0 ]; then
            #   echo "One of the vectors has a zero magnitude, skipping..."
            #    continue
            #fi

            # 计算夹角余弦值
            cos_theta=$(echo "scale=10; $dot_product / ($mod_vec1 * $mod_vec2)" | bc -l)

            # 确保夹角余弦值在有效范围内（[-1, 1]）
            cos_theta=$(echo "scale=10; if ($cos_theta < -1) -1 else if ($cos_theta > 1) 1 else $cos_theta" | bc)

            # 计算夹角（弧度制转换为角度）
            theta=$(echo "scale=10; a(${cos_theta}/sqrt(1-${cos_theta}^2)) * 360 / (4*a(1))" | bc -l)

            # 输出结果
            echo "夹角是: ${theta%.*} 度"

            # 重置捕捉状态
            found=0
        fi
    fi

done < "$file"
