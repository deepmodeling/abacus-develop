#!/bin/bash

total_time=0

# 遍历文件夹下的所有子文件夹
list=`cat CASES` 
for dir in $list; do
    # 如果 resultref.txt 文件存在
    echo dir: $dir
    if [[ -f "${dir}/result.ref" ]]; then
        # 提取含有 "time" 的行，然后提取时间值
        time_value=$(grep "time" "${dir}/result.ref" | awk '{print $2}')
        echo time_value: $time_value
        # 累加时间
        total_time=$(echo "$total_time + $time_value" | bc)
    fi
done

echo "Total time: $total_time"