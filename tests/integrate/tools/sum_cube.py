#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) < 2:
        print(f"Can't find {sys.argv[1] if len(sys.argv) > 1 else 'file'} !")
        sys.exit(1)

    input_file = sys.argv[1]

    try:
        with open(input_file, 'r') as inp:
            # skip the first two lines
            inp.readline()
            inp.readline()

            # read the 3rd line: number of atoms + origin coordinates
            line = inp.readline().split()
            if not line:
                return
            natom = int(line[0])
            # origin_x = float(line[1])
            # origin_y = float(line[2])
            # origin_z = float(line[3])

            # read the grid vectors (support non-orthogonal)
            line = inp.readline().split()
            nx = int(line[0])
            v1 = [float(line[1]), float(line[2]), float(line[3])]

            line = inp.readline().split()
            ny = int(line[0])
            v2 = [float(line[1]), float(line[2]), float(line[3])]

            line = inp.readline().split()
            nz = int(line[0])
            v3 = [float(line[1]), float(line[2]), float(line[3])]

            # calculate the volume element |v1 · (v2 × v3)|
            val0 = v2[1] * v3[2] - v2[2] * v3[1]
            val1 = v2[0] * v3[2] - v2[2] * v3[0]
            val2 = v2[0] * v3[1] - v2[1] * v3[0]
            
            volume = abs(v1[0] * val0 - v1[1] * val1 + v1[2] * val2)

            # skip the atom coordinates
            # natom can be negative in cube files sometimes?
            # C++ code: for (int i = 0; i < natom; ++i)
            # If natom is negative, loop doesn't run.
            
            atoms_to_skip = natom
            if atoms_to_skip < 0:
                 # In some cube formats, negative natom means second line of header contains units or something?
                 # Standard: "If the number of atoms is negative, that indicates that the file contains input lines... One line for each non-zero E value."
                 # But the C++ code simple loops < natom. If natom < 0, it skips 0 lines.
                 # We will mimic C++ behavior assuming it works for the files they have.
                 atoms_to_skip = 0 # Loop won't run if natom < 0

            for _ in range(atoms_to_skip):
                inp.readline()

            nr = nx * ny * nz
            
            total_sum = 0.0
            count = 0
            
            # Read grid values
            # iterate over remaining lines to handle values spread across lines
            for line in inp:
                parts = line.split()
                for part in parts:
                    total_sum += float(part)
                    count += 1
                    if count >= nr:
                        break
                if count >= nr:
                    break
            
            ne = total_sum * volume
            # cout default precision is 6, but setprecision(10) changes it.
            # Python's default float printing is usually sufficient, but let's use formatting to be sure.
            # {:.10g} prints up to 10 significant digits.
            print(f"{ne:.10g}")

    except FileNotFoundError:
        print(f"Can't find {input_file} !")
        sys.exit(1)
    except Exception as e:
        # C++ doesn't print other errors generally, but let's be safe.
        # But to be functionally equivalent, maybe we shouldn't.
        # The C++ code crashes or behaves weirdly on bad input.
        # Let's just exit 1.
        sys.exit(1)

if __name__ == "__main__":
    main()
