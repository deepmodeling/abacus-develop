Bohr2Ang = 0.529177208

with open('FORCE_CONSTANTS', 'r') as f1:
    fc = f1.readlines()
with open('FORCE_CONSTANTS_2ND', 'w') as f2:
    for idx, ii in enumerate(fc):
        char = ii.split()
        if len(char) == 3:
            for idy, jj in enumerate(char):
                char[idy] = str(float(jj) / Bohr2Ang)
            bloks = ' '
            fc[idx] = bloks.join(char)
            fc[idx] = fc[idx] + '\n'
        f2.write(fc[idx])
