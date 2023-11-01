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
        self._indices, self.dipole = self.read(self.dipolefile)
        self.dt = dt  # in fs

    @classmethod
    def read(cls, filename: PathLike):
        """Read dipole data file

        :params filename: string of dipole data file
        """
        data = np.loadtxt(filename, dtype=float)
        indices = data[:, 0]
        dipole = data[:, 1:].transpose()
        for i in range(3):
            dipole[i,:] = dipole[i,:] - dipole[i,0]

        return indices, dipole

    @property
    def dipole_data(self):
        return self.dipole

    def plot_dipole(self, fig: Figure, ax: axes.Axes, directions: list = [0, 1, 2], time_range: list = []):
        """Plot dipole data in x,y,z directions

        :params fig: (matplotlib.figure.Figure)
        :params ax: (matplotlib.axes.Axes)
        :params directions: (list) 0->X, 1->Y, 2->Z
        :params time_range: (list) range of time (in unit fs) to plot
        """

        labels = {0: 'X', 1: 'Y', 2: 'Z'}

        for direc in directions:
            ax.plot(self._indices*self.dt, self.dipole_data[direc], label=labels[direc])

        ax.set_xlabel('Times (fs)')
        ax.set_ylabel('Dipole')
        ax.legend()
        if time_range:
            ax.set_xlim(time_range)

        return fig, ax

class Efield:
    """Parse Efield Data"""

    def __init__(self, Efile: list[list[PathLike]], dt: float, nstep:int) -> None:
        self.Efile = Efile
        self.dt=dt
        self.efield=self.reade(self.Efile, nstep)

    @classmethod
    def reade(cls, Efile: list[list[PathLike]], nstep: int):
        """Read dipole data file

        :params Efile: string type 2D list of Efield data file, Efile[i][j] is the jth Efield data file in ith direction
        :params nstep: number of steps in simulation
        """
        Efield = np.zeros([3,nstep])
        for i in range(len(Efile)):
            for file in Efile[i]:
                Edata = np.loadtxt(file, dtype=float)
                Efield[i,0:(len(Edata))] += Edata[:, 1]
        return Efield
    
    @property
    def efield_data(self):
        return self.efield
    
class Absorption(Dipole, Efield):
    """Calculate Absorption Spectrum under light field"""
    def __init__(self, dipolefile: PathLike,Efile: list[list[PathLike]], dt: float) -> None:
        Dipole.__init__(self, dipolefile, dt)
        Efield.__init__(self, Efile, dt, len(self._indices))

    def padding(self, data: np.ndarray):
        """Zero padding for FFT

        :params data: (np.ndarray) data to be padded
        """
        #mask part
        Ndim=len(self._indices)*10
        index=np.linspace(0,Ndim-1,Ndim)
        t=self._indices*self.dt
        b=5
        mask=np.exp(-b*t/t[-1])
        #padding part
        padding=np.zeros(Ndim-len(data))
        data_pass=np.concatenate((data*mask, padding))
        return data_pass
    
    def alpha(self,dirc: int):
        """Calculate alpha

        :params dirc: (int) 0->X, 1->Y, 2->Z
        """
        #FFT part
        dipole=self.padding(self.dipole_data[dirc])
        efield=self.padding(self.efield_data[dirc])
        dipole_fft = np.fft.fft(dipole)
        efield_fft = np.fft.fft(efield)
        numerator=dipole_fft
        denominator=efield_fft
        alpha = np.abs((numerator/denominator).imag)
        return alpha
    
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
        Ndim=len(self._indices)*10
        index=np.linspace(0,Ndim-1,Ndim)
        energies = Freq2eV*index/self.dt/len(index)
        x_data = energies if unit == 'eV' else sc.nu2lambda(
            sc.eV/sc.h*energies)*1e9
        
        #plot the adsorption spectra and output the data
        adsorption_spectra_data = x_data[:, np.newaxis]
        for direc in directions:
            alpha = self.alpha(direc)
            ax.plot(x_data, alpha, label=labels[direc])
            adsorption_spectra_data = np.concatenate((adsorption_spectra_data, alpha[:, np.newaxis]),axis=1)
        np.savetxt('absorpation_spectra.dat', adsorption_spectra_data)

        xlabel = 'Energy (eV)' if unit == 'eV' else 'Wave Length (nm)'
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Absportion')
        ax.legend()
        if x_range:
            ax.set_xlim(x_range)
            lim_range=index[(x_data>x_range[0])&(x_data<x_range[1])]
            ax.set_ylim([0, 1.2*max(alpha[int(lim_range[0]):int(lim_range[-1])])])
        return fig, ax


if __name__ == "__main__":
    dipolefile = './SPIN1_DIPOLE'
    dipole = Dipole(dipolefile, dt=0.0024)
    Efile=[[],[],["efield_0.dat"]]
    Abs = Absorption(dipolefile, Efile, dt=0.0024)

    import matplotlib.pyplot as plt
    fig1, ax1 = plt.subplots()
    fig1, ax1 = dipole.plot_dipole(fig1, ax1)
    fig1.savefig('dipole.png')
    
    fig2, ax2 = plt.subplots()
    x_range = [0, 400]
    fig2, ax2 = Abs.plot_abs(
        fig2, ax2, x_range=x_range, unit='nm')
    fig2.savefig('abs.png')
