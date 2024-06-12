import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

current_dir = os.getcwd()

if os.path.exists(os.path.join(current_dir, "ElecStaticPot.cube")):
    input_file = os.path.join(current_dir, "ElecStaticPot.cube")
    output_file = os.path.join(current_dir, "ElecStaticPot_AVE")
    output_png = os.path.join(current_dir, "ElecStaticPot-vs-Z.png")
else:
    for dir_name in os.listdir(current_dir):
        dir_path = os.path.join(current_dir, dir_name)

        if os.path.isdir(dir_path) and dir_name.startswith("OUT."):
            input_file = os.path.join(dir_path, "ElecStaticPot.cube")
            output_file = os.path.join(dir_path, "ElecStaticPot_AVE")
            output_png = os.path.join(dir_path, "ElecStaticPot-vs-Z.png") 
            if os.path.exists(input_file):
                print(f"Processeding: {input_file}")
                break
            else:
                print(f"File does not exist: {input_file}")
                sys.exit()
                        
with open(input_file, 'r') as inpt:
    temp = inpt.readlines()

Ry_to_eV = 13.605698066819
bohr_to_Ang = 1.8897259885789

natom = int(temp[2].split()[0])
nx = int(temp[3].split()[0])
ny = int(temp[4].split()[0])
nz = int(temp[5].split()[0])
step_z_bohr = float(temp[5].split()[3])
step_z = step_z_bohr / bohr_to_Ang
z_Ang = np.arange(nz) * step_z

nrow = nz // 6
if nz % 6 != 0: nrow += 1

data_start = 6 + natom
average = np.zeros(nz)
average_eV = np.zeros(nz)

nxy = nx * ny
for ixy in range(nxy):
    zpiece = []
    for irow in range(nrow):
        zpiece += [float(_) for _ in temp[data_start + ixy * nrow + irow].split()]
    average += np.array(zpiece)

average /= nxy

average_eV = average * Ry_to_eV

with open(output_file, 'w') as output:
    output.write("Average electrostatic potential along z axis\n")
    for i in range(1, data_start):
        output.write(temp[i])
    
    output.write("iz\t\t Average(Ry)\t\t z(Angstrom)\t\t Average(eV)\n")
    for iz in range(nz):
        output.write("{:<7d}\t{:>16.9e}\t{:>16.9e}\t{:>16.9e}\n".format(iz, z_Ang[iz], average[iz], average_eV[iz]))


z_values = z_Ang
interpolation_func = interp1d(z_values, average_eV, kind='cubic')

z_interpolated = np.linspace(z_values.min(), z_values.max(), nz * 5)
average_interpolated = interpolation_func(z_interpolated)

plt.figure(figsize=(4.5, 8)) 
plt.plot(average_interpolated, z_interpolated)
plt.xlabel('Average Electrostatic Potential (eV)')
plt.ylabel('Z axis (Angstrom)')
plt.title('Average Electrostatic Potential along Z')
plt.grid(True)

plt.savefig(output_png, dpi=300)
plt.show()
