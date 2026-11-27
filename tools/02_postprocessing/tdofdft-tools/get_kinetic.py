#!/usr/bin/env python3
"""Extract the last proton's kinetic energy from MD STRU frames.

Usage: python3 get_kinetic.py [stru_file] [output_file]
Default input: OUT.autotest/STRU (out_stru=1, out_freq_ion=1, out_app_flag=1).
"""

import sys
from pathlib import Path


stru_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("OUT.autotest/STRU")
output_file = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("kinetic.txt")
md_dt = 0.0005  # fs; must match md_dt in INPUT
m_H = 1.008  # amu

# Match the constants used by ABACUS to convert STRU atomic-unit velocities.
angstrom_au = 1.8897270
au_to_fs = 2.418884326505e-2
au_velocity_to_ang_fs = 1.0 / angstrom_au / au_to_fs
amu_to_kg = 1.66053906660e-27
ang_fs_to_m_s = 1e5
eV_to_J = 1.602176634e-19
kinetic_factor = amu_to_kg * ang_fs_to_m_s**2 / eV_to_J


def proton_velocity(frame):
    data = [line.split("#", 1)[0].strip() for line in frame]
    data = [line for line in data if line]
    position_index = data.index("ATOMIC_POSITIONS") + 2
    last_velocity = None
    while position_index < len(data):
        label = data[position_index]
        count = int(data[position_index + 2])
        position_index += 3
        for atom_line in data[position_index:position_index + count]:
            fields = atom_line.split()
            if label == "H" and "v" in fields:
                velocity_index = fields.index("v")
                last_velocity = tuple(float(value) for value in fields[velocity_index + 1:velocity_index + 4])
        position_index += count
    if last_velocity is None:
        raise ValueError("Frame has no H atom velocity")
    return last_velocity


if not stru_file.is_file():
    raise SystemExit(f"Error: {stru_file} not found")

frames = []
step = None
frame = []
for line in stru_file.read_text().splitlines():
    if line.startswith("# MD step "):
        if step is not None:
            frames.append((step, frame))
        step = int(line.split()[3])
        frame = []
    elif step is not None:
        frame.append(line)
if step is not None:
    frames.append((step, frame))

if not frames:
    raise ValueError(f"No MD frames found in {stru_file}")

with output_file.open("w") as output:
    for step, frame in frames:
        velocity = proton_velocity(frame)
        velocity_squared = sum((value * au_velocity_to_ang_fs) ** 2 for value in velocity)
        time_fs = (step - 1) * md_dt
        energy = 0.5 * m_H * velocity_squared * kinetic_factor
        output.write(f"{time_fs:.4f}\t{energy:.10f}\n")

print(f"Done! {len(frames)} steps written to {output_file}")
