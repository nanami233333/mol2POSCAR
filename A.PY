import os
import glob

def convert_mol_to_poscar(input_pattern="molecule_*.mol",
                          output_dir="output",
                          margin=10.0):
    """
    将符合 input_pattern 的所有 .mol 文件转换为 POSCAR 格式，并写到 output_dir 下。
    - margin: 在 x、y、z 三个方向的原子最小/最大值区间上，各额外加的真空层 (Å 单位)。
    """
    # 定义 POSCAR 中的元素顺序（如果你的 .mol 文件还包含其他元素，请在这里补充并保证顺序统一）
    elements = ["O", "C", "H", "N", "F"]

    os.makedirs(output_dir, exist_ok=True)

    for filepath in glob.glob(input_pattern):
        # 1) 先读取 .mol 文件，将所有原子按 “元素 -> 坐标列表” 存储
        coords_by_elem = {el: [] for el in elements}
        # 同时记录所有坐标以便后面计算最小/最大值
        all_x, all_y, all_z = [], [], []

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Molecule 文件通常前 4 行是 header，从第 5 行（索引 4）开始才是真正的原子坐标
        # 遇到以 'M' 开头的行（如 "M  END"）就停止
        for line in lines[4:]:
            stripped = line.strip()
            if stripped.startswith("M"):
                break

            parts = stripped.split()
            if len(parts) >= 4:
                x_f = float(parts[0])
                y_f = float(parts[1])
                z_f = float(parts[2])
                el = parts[3]
                # 只收集属于我们 elements 列表里的元素，其他类型忽略
                if el in coords_by_elem:
                    coords_by_elem[el].append((x_f, y_f, z_f))
                    all_x.append(x_f)
                    all_y.append(y_f)
                    all_z.append(z_f)
                else:
                    # 如果你的 .mol 里还有别的元素，请把它加到 elements 里；
                    # 这里先简单提示，然后跳过
                    print(f"警告：在文件 {filepath} 中发现未列入 {elements} 的元素 '{el}'，已忽略。")

        # 如果没有任何原子读到，则跳过
        if len(all_x) == 0:
            print(f"注意：文件 {filepath} 中没有读取到任何支持的原子，已跳过。")
            continue

        # 2) 计算 xyz 三个方向的最小/最大值
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        z_min, z_max = min(all_z), max(all_z)

        # 3) 根据 (max - min) 加上 margin，得到正交晶胞的边长
        Lx = (x_max - x_min) + margin
        Ly = (y_max - y_min) + margin
        Lz = (z_max - z_min) + margin

        # 4) 计算每个原子新的 “居中 + 平移” 后的坐标
        #    公式： new_x = x - (x_min + x_max)/2 + Lx/2 = x - x_min + margin/2
        #    下面遍历每个元素列表，把坐标转换好再一起写入 POSCAR。
        shifted_by_elem = {el: [] for el in elements}
        for el in elements:
            for (x_f, y_f, z_f) in coords_by_elem[el]:
                x_new = x_f - x_min + margin / 2.0
                y_new = y_f - y_min + margin / 2.0
                z_new = z_f - z_min + margin / 2.0
                shifted_by_elem[el].append((x_new, y_new, z_new))

        # 5) 准备 POSCAR 文件内容
        poscar_lines = []
        # (1) 第一行：可以留空，也可以写注释。例如：# molecule_123.mol 
        poscar_lines.append("")  

        # (2) 第二行：缩放因子，一般写 1.0
        poscar_lines.append("1.00000000000000")

        # (3)(4)(5) 晶格向量（正交晶胞）
        poscar_lines.append(f"{Lx:.6f} 0.000000 0.000000")
        poscar_lines.append(f"0.000000 {Ly:.6f} 0.000000")
        poscar_lines.append(f"0.000000 0.000000 {Lz:.6f}")

        # (6) 元素顺序
        poscar_lines.append(" ".join(elements))

        # (7) 各元素的计数
        counts = [len(coords_by_elem[el]) for el in elements]
        poscar_lines.append(" ".join(str(c) for c in counts))

        # (8) 坐标类型：Cartesian
        poscar_lines.append("Cartesian")

        # (9) 坐标部分：一定要“逐元素顺序写入”，保证和上面第 (6)(7) 行的顺序一致
        for el in elements:
            for (x_new, y_new, z_new) in shifted_by_elem[el]:
                poscar_lines.append(f"{x_new:.6f} {y_new:.6f} {z_new:.6f} {el}")

        # 6) 写入到对应的 .POSCAR 文件里
        filename = os.path.basename(filepath)     # e.g. "molecule_123.mol"
        name, _ = os.path.splitext(filename)      # e.g. "molecule_123"
        out_path = os.path.join(output_dir, f"{name}.POSCAR")
        with open(out_path, "w") as fout:
            fout.write("\n".join(poscar_lines))

        print(f"生成成功：{out_path}")

if __name__ == "__main__":
    # margin 可根据实际需要增大/缩小，比如改成 8.0 或 12.0：
    convert_mol_to_poscar(input_pattern="molecule_*.mol", output_dir="output", margin=10.0)
