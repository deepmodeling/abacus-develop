def read_mmn_data(filename):
    """读取 .mmn 文件，跳过 header，返回数据行列表（每行两个 float）"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    # 跳过前3行（通常是 header）
    data_lines = lines[3:]
    data = []
    for i, line in enumerate(data_lines, start=4):  # 原始文件行号从4开始
        parts = line.split()
        if len(parts) >= 2:
            try:
                real = float(parts[0])
                imag = float(parts[1])
                data.append((real, imag, i))  # 保存值和原始行号
            except ValueError:
                # 跳过无法解析的行（如空行）
                continue
    return data

def main():
    file1 = 'OUT.autotest/diamond.mmn'
    file2 = 'OUT.develop/diamond.mmn'

    data1 = read_mmn_data(file1)
    data2 = read_mmn_data(file2)

    if len(data1) != len(data2):
        print("警告：两个文件数据行数不一致！")
        min_len = min(len(data1), len(data2))
    else:
        min_len = len(data1)

    max_diff = -1.0
    max_info = None  # (diff, file1_val, file2_val, line_num, col)

    for i in range(min_len):
        r1, i1, line1 = data1[i]
        r2, i2, line2 = data2[i]
        # 实部差值
        dr = abs(r2 - r1)
        # 虚部差值
        di = abs(i2 - i1)

        if dr > max_diff:
            max_diff = dr
            max_info = (dr, r1, r2, line1, 'real')
        if di > max_diff:
            max_diff = di
            max_info = (di, i1, i2, line1, 'imag')

    if max_info:
        diff_val, val1, val2, line_num, part = max_info
        print("🔍 最大绝对偏差详情:")
        print(f"  文件1值: {val1:.12e}")
        print(f"  文件2值: {val2:.12e}")
        print(f"  差值（绝对值）: {diff_val:.3e}")
        print(f"  所在行号（在原始 .mmn 文件中）: {line_num}")
        print(f"  分量: {part} 部分")
        print(f"\n📌 建议检查 {file1} 和 {file2} 的第 {line_num} 行")
    else:
        print("未找到有效数据。")

    # 可选：输出所有差值到文件
    with open('mmn_full_diff.txt', 'w') as out:
        out.write("# line  real1  real2  real_diff  imag1  imag2  imag_diff\n")
        for i in range(min_len):
            r1, i1, ln = data1[i]
            r2, i2, _ = data2[i]
            dr = r2 - r1
            di = i2 - i1
            out.write(f"{ln}  {r1:.12e}  {r2:.12e}  {dr:.3e}  {i1:.12e}  {i2:.12e}  {di:.3e}\n")
    print("\n📊 完整差值已保存到 'mmn_full_diff.txt'")
    
if __name__ == '__main__':
    main()