import numpy as np

def check_if_reduce(path = "./"):

    if path[-1] != "/":
        path += "/"
    log_fname = path + f"kpoints"
    line = "start"

    with open(log_fname, "r") as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip()
        if line.startswith("KPT     DIRECT_X"):
            return True
    return False

def parse(path = "./"):
    if check_if_reduce(path):
        return parse_symmetry_true(path)
    else:
        return parse_symmetry_false(path)

def parse_symmetry_true(path = "./"):

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

def parse_symmetry_false(path = "./"):

    if path[-1] != "/":
        path += "/"
    log_fname = path + f"kpoints"
    line = "start"

    kpoints = []
    nkpts = 0
    with open(log_fname, "r") as f:
        contents = f.readlines()
    for line in contents:
        line = line.strip()
        if line.startswith("nkstot now = "):
            nkpts = int(line.split()[-1])
        if line[0].isdigit():
            kpoint = np.array([float(x) for x in line.split()[1:4]])
            kpoints.append(kpoint)
            nkpts -= 1
        if nkpts == 0:
            break
    return kpoints, kpoints

if __name__ == "__main__":
    kpoints, equiv_kpoints = parse("./work")
    print(kpoints)
    print(equiv_kpoints)