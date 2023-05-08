from abacus_plot.dipole import Dipole
from abacus_plot.dipole import Absorption
import matplotlib.pyplot as plt

dipolefile = './SPIN1_DIPOLE'
dipole = Dipole(dipolefile, dt=0.0024)
Efile=[[],[],["efield_0.dat","efield_1.dat"]]
Abs = Absorption(dipolefile, Efile, dt=0.0024)

fig1, ax1 = plt.subplots()
fig1, ax1 = dipole.plot_dipole(fig1, ax1)
fig1.savefig('dipole.png')

fig2, ax2 = plt.subplots()
x_range = [0, 20]
fig2, ax2 = Abs.plot_abs(
fig2, ax2, x_range=x_range,directions=[2], unit='eV')
fig2.savefig('abs.png')