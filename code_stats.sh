#!/bin/bash
#统计.cpp文件
cpp_count=$(find source -name "*.cpp" | wc -l)
cpp_lines=$(find source -name "*.cpp" | xargs cat 2>/dev/null | wc -l)
cpp_zhu=$(find source -name "*.cpp" | xargs cat 2>/dev/null | grep -E "^[[:space:]]*(//|/\*|\*|.*\*/)" | wc -l)
#统计.h文件
h_count=$(find source -name "*.h" | wc -l)
h_lines=$(find source -name "*.h" | xargs cat 2>/dev/null | wc -l)
h_zhu=$(find source -name "*.h" | xargs cat 2>/dev/null | grep -E "^[[:space:]]*(//|/\*|\*|.*\*/)" | wc -l)
#分别计算注释率
cpprate=$(echo "scale=2; 100 *  $cpp_zhu / $cpp_lines " | bc)
hrate=$(echo "scale=2; 100 *  $h_zhu / $h_lines " | bc)
echo ".cpp 文件数量: $cpp_count"
echo ".cpp 总行数: $cpp_lines"
echo ".cpp 注释行数: $cpp_zhu"
echo ".cpp 注释率：${cpprate}%"
echo ".h 文件数量: $h_count"
echo ".h 总行数: $h_lines"
echo ".h 注释行数: $h_zhu"
echo ".h 注释率：${hrate}%"

