import numpy as np
from os import PathLike
from matplotlib.figure import Figure
from matplotlib import axes
import scipy.constants as sc

Freq2eV = sc.h/sc.eV*1e15  # 1/fs to eV


class Dipole:
    """Parse Dipole Data"""

    def __init__(self, dipolefile: PathLike, dt: float) -> None:
        self.dipolefile = dipolefile
        self._indices, self.dipole_x, self.dipole_y, self.dipole_z = self.read(
            self.dipolefile)
        self.dt = dt  # in fs
        self.time = self._indices*dt
        self.energies = Freq2eV*self._indices/dt/len(self._indices)

    @classmethod
    def read(cls, filename: PathLike):
        """Read dipole data file

        :params filename: string of dipole data file
        """

        data = np.loadtxt(filename, dtype=float)
        indices = data[:, 0]
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]

        return indices, x, y, z

    @property
    def dipole_data(self):
        return {0: self.dipole_x, 1: self.dipole_y, 2: self.dipole_z}

    def plot_dipole(self, fig: Figure, ax: axes.Axes, directions: list = [0, 1, 2], time_range: list = []):
        """Plot dipole data in x,y,z directions

        :params fig: (matplotlib.figure.Figure)
        :params ax: (matplotlib.axes.Axes)
        :params directions: (list) 0->X, 1->Y, 2->Z
        :params time_range: (list) range of time (in unit fs) to plot
        """

        labels = {0: 'X', 1: 'Y', 2: 'Z'}

        for direc in directions:
            ax.plot(self.time, self.dipole_data[direc], label=labels[direc])

        ax.set_xlabel('Times (fs)')
        ax.set_ylabel('Dipole')
        ax.legend()
        if time_range:
            ax.set_xlim(time_range)

        return fig, ax

    @property
    def alpha_x(self):
        return np.fft.fft(self.dipole_x)

    @property
    def alpha_y(self):
        return np.fft.fft(self.dipole_y)

    @property
    def alpha_z(self):
        return np.fft.fft(self.dipole_z)

    @property
    def alpha_data(self):
        return {0: self.alpha_x, 1: self.alpha_y, 2: self.alpha_z}

    def get_abs(self, direc: int):
        S = np.abs(self.alpha_data[direc].imag)
        return S

    def plot_abs(self, fig: Figure, ax: axes.Axes, directions: list = [0, 1, 2], x_range: list = [], unit: str = 'eV'):
        """Plot Absportion Spectrum under Delta light field in x,y,z directions

        :params fig: (matplotlib.figure.Figure)
        :params ax: (matplotlib.axes.Axes)
        :params directions: (list) 0->X, 1->Y, 2->Z
        :params x_range: (list) range of energies (in unit eV) to plot
        :params unit: (str) 
        """

        assert unit in ['eV', 'nm']
        labels = {0: 'X', 1: 'Y', 2: 'Z'}
        x_data = self.energies if unit == 'eV' else sc.nu2lambda(
            sc.eV/sc.h*self.energies)*1e9
        
        #plot the adsorption spectra and output the data
        adsorption_spectra_data = x_data[:, np.newaxis]
        for direc in directions:
            ax.plot(x_data, self.get_abs(direc), label=labels[direc])
            adsorption_spectra_data = np.concatenate((adsorption_spectra_data, self.get_abs(direc)[:, np.newaxis]),axis=1)
            if direc != directions[-1]:
                adsorption_spectra_data = np.concatenate((adsorption_spectra_data, x_data[:, np.newaxis]),axis=1)
        np.savetxt('absorpation_spectra.dat', adsorption_spectra_data)

        xlabel = 'Energy (eV)' if unit == 'eV' else 'Wave Length (nm)'
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Absportion')
        ax.legend()
        if x_range:
            ax.set_xlim(x_range)

        return fig, ax


if __name__ == "__main__":
    dipolefile = './SPIN1_DIPOLE'
    dipole = Dipole(dipolefile, dt=0.0034)

    import matplotlib.pyplot as plt
    fig1, ax1 = plt.subplots()
    fig1, ax1 = dipole.plot_dipole(fig1, ax1)
    fig1.savefig('dipole.png')

    fig2, ax2 = plt.subplots()
    x_range = [0, 400]
    fig2, ax2 = dipole.plot_abs(
        fig2, ax2, x_range=x_range, unit='nm')
    fig2.savefig('abs.png')
