def read_amn_data(filename):
    """读取 .amn 文件，跳过 header，返回数据行列表"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    # 跳过前2行（header）
    data_lines = lines[2:]
    data = []
    for i, line in enumerate(data_lines, start=3):  # 原始文件行号从3开始
        parts = line.split()
        if len(parts) >= 5:
            try:
                band_idx = int(parts[0])
                wannier_idx = int(parts[1])
                kpt_idx = int(parts[2])
                real = float(parts[3])
                imag = float(parts[4])
                data.append((band_idx, wannier_idx, kpt_idx, real, imag, i))
            except (ValueError, IndexError):
                continue
    return data

def main():
    file1 = 'OUT.autotest/diamond.amn'
    file2 = 'OUT.develop/diamond.amn'

    data1 = read_amn_data(file1)
    data2 = read_amn_data(file2)

    if len(data1) != len(data2):
        print(f"警告：两个文件数据行数不一致！{len(data1)} vs {len(data2)}")
        min_len = min(len(data1), len(data2))
    else:
        min_len = len(data1)

    max_diff = -1.0
    max_info = None  # (diff, val1, val2, line, part, indices)

    for i in range(min_len):
        b1, w1, k1, r1, i1, ln1 = data1[i]
        b2, w2, k2, r2, i2, ln2 = data2[i]

        # 验证索引是否一致（应一致）
        if (b1, w1, k1) != (b2, w2, k2):
            print(f"警告：第 {ln1} 行索引不匹配: ({b1},{w1},{k1}) vs ({b2},{w2},{k2})")
            continue

        dr = abs(r2 - r1)
        di = abs(i2 - i1)

        if dr > max_diff:
            max_diff = dr
            max_info = (dr, r1, r2, ln1, 'real', (b1, w1, k1))
        if di > max_diff:
            max_diff = di
            max_info = (di, i1, i2, ln1, 'imag', (b1, w1, k1))

    if max_info:
        diff_val, val1, val2, line_num, part, (b, w, k) = max_info
        print("🔍 最大绝对偏差详情:")
        print(f"  文件1值: {val1:.12e}")
        print(f"  文件2值: {val2:.12e}")
        print(f"  差值（绝对值）: {diff_val:.3e}")
        print(f"  所在行号（在原始 .amn 文件中）: {line_num}")
        print(f"  分量: {part} 部分")
        print(f"  索引: band={b}, wannier={w}, kpoint={k}")
        print(f"\n📌 建议检查 {file1} 和 {file2} 的第 {line_num} 行")
    else:
        print("未找到有效数据。")

    # 可选：输出完整差值
    with open('amn_full_diff.txt', 'w') as out:
        out.write("# line band wannier kpt real1 real2 real_diff imag1 imag2 imag_diff\n")
        for i in range(min_len):
            b, w, k, r1, i1, ln = data1[i]
            _, _, _, r2, i2, _ = data2[i]
            dr = r2 - r1
            di = i2 - i1
            out.write(f"{ln} {b} {w} {k} {r1:.6e} {r2:.6e} {dr:.3e} {i1:.6e} {i2:.6e} {di:.3e}\n")
    print("\n📊 完整差值已保存到 'amn_full_diff.txt'")
    
if __name__ == '__main__':
    main()