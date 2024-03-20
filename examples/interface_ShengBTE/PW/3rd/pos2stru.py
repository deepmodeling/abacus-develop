import glob
from ase.io import read, write
from pathlib import Path

pos = glob.glob('3RD.POSCAR.*')
for ii in pos:
    idx = ii.split('.')[-1]
    print(ii, idx)
    cs_dir = './'
    cs_vasp = Path(cs_dir, ii)
    cs_atoms = read(cs_vasp, format='vasp')
    stru = 'STRU_' + str(idx)
    print(stru)
    cs_stru = Path(cs_dir, stru)
    pp = {'Si':'Si_ONCV_PBE-1.2.upf'}
    write(cs_stru, cs_atoms, format='abacus', pp=pp)
    # basis = {'Si':'Si_gga_7au_100Ry_2s2p1d.orb'}
    # write(cs_stru, cs_atoms, format='abacus', pp=pp, basis=basis)
    