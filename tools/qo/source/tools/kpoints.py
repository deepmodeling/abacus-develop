import numpy as np

def parse(path = "./"):

    if path[-1] != "/":
        path += "/"
    log_fname = path + f"kpoints"
    line = "start"

    kpoints = []
    equiv_kpoints = []

    read_kpoint_reduction_information = False
    with open(log_fname, "r") as f:
        while len(line) > 0:
            line = f.readline().strip()
            
            if line.startswith("KPT     DIRECT_X"):
                read_kpoint_reduction_information = True
                continue
            if line.startswith("nkstot = "):
                nkpts = line.split()[2]

            if read_kpoint_reduction_information:
                
                kpt_ind, x, y, z, irr_kpt_ind, _x, _y, _z = line.split()
                if int(irr_kpt_ind) > len(kpoints):
                    kpoint = np.array([float(_x), float(_y), float(_z)])
                    kpoints.append(kpoint)
                    equiv_kpoints.append([kpoint])
                else:
                    kpoint = np.array([float(x), float(y), float(z)])
                    equiv_kpoints[int(irr_kpt_ind)-1].append(kpoint)
                
                if int(kpt_ind) == int(nkpts):
                    break
                continue


    return kpoints, equiv_kpoints

if __name__ == "__main__":
    kpoints, equiv_kpoints = parse("./work")
    print(kpoints)
    print(equiv_kpoints)