import numpy as np

input_file = "ElecStaticPot.cube"
output_file = "ElecStaticPot_AVE"

with open(input_file, 'r') as inpt:
    temp = inpt.readlines()

natom = int(temp[2].split()[0])
nx = int(temp[3].split()[0])
ny = int(temp[4].split()[0])
nz = int(temp[5].split()[0])

nrow = nz // 6
if nz % 6 != 0: nrow += 1

data_start = 6 + natom
average = np.zeros(nz)

nxy = nx * ny
for ixy in range(nxy):
    zpiece = []
    for irow in range(nrow):
        zpiece += [float(_) for _ in temp[data_start + ixy * nrow + irow].split()]
    average += np.array(zpiece)

average /= nxy

with open(output_file, 'w') as output:
    output.write("Average electrostatic potential along z axis\n")
    for i in range(1, data_start):
        output.write(temp[i])
    
    output.write("iz\t\taverage\n")
    for iz in range(nz):
        output.write("{0}\t\t{1:.9e}\n".format(iz, average[iz]))
