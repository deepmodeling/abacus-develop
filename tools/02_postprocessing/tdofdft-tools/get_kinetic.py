#!/usr/bin/env python3
"""
Extract proton kinetic energy from TDOFDFT MD_dump output.

Usage:
    python3 get_kinetic.py [md_dump_file] [output_file]

Default:
    md_dump_file = OUT.autotest/MD_dump
    output_file  = kinetic.txt

Output format:
    time(fs)    kinetic_energy(eV)
"""

import sys
import os

# ========== Parameters ==========
md_dump_file = "OUT.autotest/MD_dump"
output_file = "kinetic.txt"
md_dt = 0.0005          # fs, must match md_dt in INPUT
m_H = 1.008             # amu, proton mass
# ================================

if len(sys.argv) >= 2:
    md_dump_file = sys.argv[1]
if len(sys.argv) >= 3:
    output_file = sys.argv[2]

# Conversion factor: E_k(eV) = 0.5 * m(amu) * v(Ang/fs)^2 * factor
amu_to_kg = 1.66053906660e-27
ang_fs_to_m_s = 1e5         # 1 Ang/fs = 10^5 m/s
eV_to_J = 1.602176634e-19
factor = amu_to_kg * ang_fs_to_m_s**2 / eV_to_J   # ≈ 103.6427

if not os.path.exists(md_dump_file):
    print(f"Error: {md_dump_file} not found!")
    sys.exit(1)

with open(md_dump_file, "r") as f:
    lines = f.readlines()

kinetic_data = []
i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith("MDSTEP:"):
        step = int(line.split(":")[1].strip())
        time_fs = step * md_dt

        i += 1
        # Skip lattice + virial + header lines to reach atom data
        while i < len(lines) and "INDEX" not in lines[i] and "MDSTEP:" not in lines[i]:
            i += 1
        i += 1  # skip INDEX header

        # Atom data: INDEX LABEL x y z Fx Fy Fz Vx Vy Vz
        # Find the last H atom (proton)
        last_vx = last_vy = last_vz = 0.0
        while i < len(lines) and not lines[i].strip().startswith("MDSTEP:"):
            parts = lines[i].split()
            if len(parts) >= 11 and parts[1] == "H":
                last_vx = float(parts[8])
                last_vy = float(parts[9])
                last_vz = float(parts[10])
            i += 1

        v_sq = last_vx**2 + last_vy**2 + last_vz**2
        ek_eV = 0.5 * m_H * v_sq * factor
        kinetic_data.append((time_fs, ek_eV))
        continue
    i += 1

with open(output_file, "w") as f:
    for t, ek in kinetic_data:
        f.write(f"{t:.4f}\t{ek:.10f}\n")

print(f"Done! {len(kinetic_data)} steps written to {output_file}")
