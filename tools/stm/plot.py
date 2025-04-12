from stm import STM
import matplotlib.pyplot as plt
import numpy as np



dirname='OUT.autotest'
bias = 2
z = 8.0

stm = STM(dirname)


# 1. Constant current 2-d scan
c = stm.get_averaged_current(bias, z)
x, y, h = stm.scan(bias, c, repeat=(3, 5))

plt.gca().axis('equal')
plt.contourf(x, y, h, 40)
plt.colorbar()
plt.savefig('2d_current.png')


# 2. Constant height 2-d scan
plt.figure()
plt.gca().axis('equal')
x, y, I = stm.scan2(bias, z, repeat=(3, 5))
plt.contourf(x, y, I, 40)
plt.colorbar()
plt.savefig('2d_height.png')


# 3. Constant current line scan
plt.figure()
a = stm.atoms.cell[0, 0]
x, y = stm.linescan(bias, c, [0, 0], [2 * a, 0])
plt.plot(x, y)
plt.savefig('line.png')


# Scanning tunneling spectroscopy
# 4. dI/dV spectrum
plt.figure()
biasstart = 2.0
biasend = 3.0
biasstep = 0.1
bias, I, dIdV = stm.sts(0, 0, z, biasstart, biasend, biasstep)
plt.plot(bias, I, label='I')
plt.plot(bias, dIdV, label='dIdV')
plt.xlim(biasstart, biasend)
plt.legend()
plt.savefig('dIdV.png')


# 5. dI/dV map
biasstart = -2.0
biasend = 2.0
biasstep = 0.1
start_points = [0, 0, 8]
end_points = [2, 0, 8]
distance = np.linalg.norm(np.array(end_points) - np.array(start_points))
biases, points, dIdV_map = stm.line_sts(biasstart, biasend, biasstep, start_points, end_points, 50)

plt.figure(figsize=(10, 6), dpi=150)

# 2D Colormap
im = plt.pcolormesh(biases, points, dIdV_map, cmap='plasma', shading='auto')
# optional cmap: 'RdYlBu', 'plasma', 'inferno', 'viridis'

# color bar
cbar = plt.colorbar(im, label='dI/dV', pad=0.02)

# labels
plt.xlabel('Sample Bias (V)', fontsize=12)
plt.ylabel('distance to the start point (Angstrom)', fontsize=12)
plt.xlim(biasstart, biasend)
plt.ylim(0, distance)
plt.title('dI/dV (~density of states) map', fontsize=14)

plt.tight_layout()
plt.savefig('dIdV_map.png', bbox_inches='tight', transparent=False)
plt.show()