import numpy as np
from source.tools.basic_functions import make_complex

def parse(nkpts: int, path = "./"):
    """read QO overlap matrix S(k) from file

    Args:
        nkpts (int): number of kpoints
    """
    qo_ovlp = []
    kpoints = []

    if path[-1] != "/":
        path += "/"
    for ik in range(nkpts):
        qo_fname = path + f"QO_ovlp_{ik}.dat"
        qo_ovlp_k = []
        with open(qo_fname, "r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("KPOINT_COORDINATE:"):
                _words = line.split(":")
                _words = _words[-1].split()
                kpoints.append([float(number) for number in _words])
            else:
                qo_ovlp_k.append([make_complex(number) for number in line.split()])
        qo_ovlp.append(np.array(qo_ovlp_k))
    
    _qo_ovlp_R0 = np.zeros_like(qo_ovlp[0])
    for ik in range(nkpts):
        _qo_ovlp_R0 += qo_ovlp[ik]
    print(_qo_ovlp_R0)
    return np.array(qo_ovlp), np.array(kpoints)

if __name__ == "__main__":
    qo_ovlp = parse(2)
    print(qo_ovlp)
    print(qo_ovlp.shape)