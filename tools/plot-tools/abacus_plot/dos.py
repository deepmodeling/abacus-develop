'''
Date: 2021-12-29 15:02:13
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:18:21
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
import numpy as np
from os import PathLike
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Union, Any
from matplotlib.figure import Figure
from matplotlib import axes

from abacus_plot.utils import (energy_minus_efermi, get_angular_momentum_label,
                               get_angular_momentum_name, remove_empty, handle_data, parse_projected_data)


class DOS:
    """Parse DOS data"""

    def __init__(self) -> None:
        self.nspin = 1
        if self.nspin in [1, 4]:
            self._nsplit = 1
        elif self.nspin == 2:
            self._nsplit = 2

    def _plot(self, dosplot, energy_f, dos, label):
        if self._nsplit == 1:
            dosplot.ax.plot(energy_f, dos,

                            lw=dosplot._lw, label=label)
        elif self._nsplit == 2:
            dos_up, dos_dw = np.split(dos, self._nsplit, axis=1)
            dosplot.ax.plot(energy_f, dos_up, lw=dosplot._lw, linestyle="-",
                            label=f'{label}'+r'$\uparrow$')
            dos_dw = -dos_dw
            dosplot.ax.plot(energy_f, dos_dw, lw=dosplot._lw, linestyle="--",
                            label=f'{label}'+r'$\downarrow$')

        return dosplot.ax

    @classmethod
    def set_vcband(cls, energy: Sequence, dos: Sequence, prec=0.01) -> Tuple[namedtuple, namedtuple]:
        """Separate valence and conduct band

        :params energy: band energy after subtracting the Fermi level
        :params dos: density of state
        :params prec: dos below this value thought to be zero. Default: 0.01
        """

        Band = namedtuple('Band', ['band', 'value'])

        # valance band
        band_vbm_index = np.where(energy <= 0)[0]
        evbm = energy[np.where(dos[band_vbm_index] > prec)[0][-1]][0]
        band_vbm = energy[band_vbm_index]
        vb = Band(band_vbm-evbm, evbm-evbm)

        # conduct band
        band_cbm_index = np.where(energy > 0)[0]
        band_cbm = energy[band_cbm_index]
        ecbm = energy[np.where(dos[band_cbm_index] > prec)[
            0][0]+band_cbm_index[0]][0]
        cb = Band(band_cbm-evbm, ecbm-evbm)

        return vb, cb

    @classmethod
    def info(cls, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        gap = cls.bandgap(vb, cb)
        print("--------------------------Density of State--------------------------", flush=True)
        print(f"{'Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap


class DOSPlot:
    """Plot density of state(DOS)"""

    def __init__(self, fig: Figure = None, ax: axes.Axes = None, **kwargs) -> None:
        self.fig = fig
        self.ax = ax
        self._lw = kwargs.pop('lw', 2)
        self._bwidth = kwargs.pop('bwdith', 3)
        self.plot_params = kwargs

    def _set_figure(self, energy_range: Sequence = [], dos_range: Sequence = [], notes: Dict = {}):
        """set figure and axes for plotting

        :params energy_range: range of energy
        :params dos_range: range of dos
        """

        # y-axis
        if dos_range:
            self.ax.set_ylim(dos_range[0], dos_range[1])
        if "ylabel_params" in self.plot_params.keys():
            self.ax.set_ylabel("DOS", **self.plot_params["ylabel_params"])
        else:
            self.ax.set_ylabel("DOS", size=25)

        # x-axis
        if energy_range:
            self.ax.set_xlim(energy_range[0], energy_range[1])

        # notes
        if notes:
            from matplotlib.offsetbox import AnchoredText
            if "s" in notes.keys() and len(notes.keys()) == 1:
                self.ax.add_artist(AnchoredText(notes["s"], loc='upper left', prop=dict(size=25),
                                                borderpad=0.2, frameon=False))
            else:
                self.ax.add_artist(AnchoredText(**self.plot_params["notes"]))

        # ticks
        if "tick_params" in self.plot_params.keys():
            self.ax.tick_params(**self.plot_params["tick_params"])
        else:
            self.ax.tick_params(labelsize=25)

        # frame
        bwidth = self._bwidth
        self.ax.spines['top'].set_linewidth(bwidth)
        self.ax.spines['right'].set_linewidth(bwidth)
        self.ax.spines['left'].set_linewidth(bwidth)
        self.ax.spines['bottom'].set_linewidth(bwidth)

        if "vline_params" in self.plot_params.keys():
            self.ax.axvline(0, **self.plot_params["vline_params"])
        else:
            self.ax.axvline(0, linestyle="--", c='b', lw=1.0)

        if "legend_prop" in self.plot_params.keys():
            self.ax.legend(prop=self.plot_params["legend_prop"])
        else:
            self.ax.legend(prop={'size': 15})


class TDOS(DOS):
    """Parse total DOS data"""

    def __init__(self, tdosfile: PathLike = None) -> None:
        super().__init__()
        self.tdosfile = tdosfile
        self._read()

    def _read(self) -> tuple:
        """Read total DOS data file

        :params tdosfile: string of TDOS data file
        """

        data = np.loadtxt(self.tdosfile)
        self.nspin = data.shape[1]-1
        if self.nspin == 1:
            self.energy, self.dos = np.split(data, self.nspin+1, axis=1)
        elif self.nspin == 2:
            self.energy, dos_up, dos_dw = np.split(data, self.nspin+1, axis=1)
            self.dos = np.hstack(dos_up, dos_dw)

    def _shift_energy(self, efermi: float = 0, shift: bool = False, prec: float = 0.01):
        if shift:
            vb, cb = self.set_vcband(
                energy_minus_efermi(self.energy, efermi), self.dos, prec)
            self.info(vb, cb)
            energy_f = np.concatenate((vb.band, cb.band))
        else:
            energy_f = energy_minus_efermi(self.energy, efermi)

        return energy_f

    def plot(self, fig: Figure, ax: Union[axes.Axes, Sequence[axes.Axes]], efermi: float = 0, energy_range: Sequence[float] = [], dos_range: Sequence[float] = [], shift: bool = False, prec: float = 0.01, **kwargs):

        energy_f = self._shift_energy(efermi, shift, prec)

        dosplot = DOSPlot(fig, ax, **kwargs)
        dosplot.ax = self._plot(dosplot, energy_f, self.dos, "TDOS")
        if "notes" in dosplot.plot_params.keys():
            dosplot._set_figure(energy_range, dos_range,
                                notes=dosplot.plot_params["notes"])
        else:
            dosplot._set_figure(energy_range, dos_range)

        return dosplot


class PDOS(DOS):
    """Parse partial DOS data"""

    def __init__(self, pdosfile: PathLike = None) -> None:
        super().__init__()
        self.pdosfile = pdosfile
        self._read()

    def _read(self):
        """Read partial DOS data file

        :params pdosfile: string of PDOS data file
        """

        from lxml import etree
        pdosdata = etree.parse(self.pdosfile)
        root = pdosdata.getroot()
        self.nspin = int(root.xpath('//nspin')[0].text.replace(' ', ''))
        norbitals = int(root.xpath('//norbitals')[0].text.replace(' ', ''))
        self.eunit = root.xpath('//energy_values/@units')[0].replace(' ', '')
        e_list = root.xpath(
            '//energy_values')[0].text.replace(' ', '').split('\n')
        remove_empty(e_list)
        self.orbitals = []
        for i in range(norbitals):
            orb = OrderedDict()
            orb['index'] = int(root.xpath(
                '//orbital/@index')[i].replace(' ', ''))
            orb['atom_index'] = int(root.xpath(
                '//orbital/@atom_index')[i].replace(' ', ''))
            orb['species'] = root.xpath(
                '//orbital/@species')[i].replace(' ', '')
            orb['l'] = int(root.xpath('//orbital/@l')[i].replace(' ', ''))
            orb['m'] = int(root.xpath('//orbital/@m')[i].replace(' ', ''))
            orb['z'] = int(root.xpath('//orbital/@z')[i].replace(' ', ''))
            data = root.xpath('//data')[i].text.split('\n')
            data = handle_data(data)
            remove_empty(data)
            orb['data'] = np.asarray(data, dtype=float)
            self.orbitals.append(orb)

        self.energy = np.reshape(e_list, newshape=(-1, 1)).astype(float)

    def _all_sum(self) -> Tuple[np.ndarray, int]:
        res = np.zeros_like(self.orbitals[0]["data"], dtype=float)
        for orb in self.orbitals:
            res = res + orb['data']
        return res

    def _write(self, species: Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], keyname='', outdir: PathLike = './'):
        """Write parsed partial dos data to files

        Args:
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            keyname (str): the keyword that extracts the PDOS. Allowed values: 'index', 'atom_index', 'species'
        """

        dos, totnum = parse_projected_data(self.orbitals, species, keyname)
        fmt = ['%13.7f', '%15.8f'] if self._nsplit == 1 else [
            '%13.7f', '%15.8f', '%15.8f']
        file_dir = Path(f"{outdir}", "PDOS_FILE")
        file_dir.mkdir(exist_ok=True)

        if isinstance(species, (list, tuple)):
            for elem in dos.keys():
                header_list = ['']
                data = np.hstack((self.energy.reshape(-1, 1), dos[elem]))
                with open(file_dir/f"{keyname}-{elem}.dat", 'w') as f:
                    header_list.append(
                        f"Partial DOS for {keyname}: {elem}")
                    header_list.append('')
                    for orb in self.orbitals:
                        if orb[keyname] == elem:
                            header_list.append(
                                f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                    header_list.append('')
                    header_list.append('\tEnergy'+10*' ' +
                                       'spin 1'+8*' '+'spin 2')
                    header_list.append('')
                    header = '\n'.join(header_list)
                    np.savetxt(f, data, fmt=fmt, header=header)

        elif isinstance(species, dict):
            for elem in dos.keys():
                elem_file_dir = file_dir/f"{keyname}-{elem}"
                elem_file_dir.mkdir(exist_ok=True)
                for ang in dos[elem].keys():
                    l_index = int(ang)
                    if isinstance(dos[elem][ang], dict):
                        for mag in dos[elem][ang].keys():
                            header_list = ['']
                            data = np.hstack(
                                (self.energy.reshape(-1, 1), dos[elem][ang][mag]))
                            m_index = int(mag)
                            with open(elem_file_dir/f"{keyname}-{elem}_{ang}_{mag}.dat", 'w') as f:
                                header_list.append(
                                    f"Partial DOS for {keyname}: {elem}")
                                header_list.append('')
                                for orb in self.orbitals:
                                    if orb[keyname] == elem and orb["l"] == l_index and orb["m"] == m_index:
                                        header_list.append(
                                            f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                                header_list.append('')
                                header_list.append(
                                    '\tEnergy'+10*' '+'spin 1'+8*' '+'spin 2')
                                header_list.append('')
                                header = '\n'.join(header_list)
                                np.savetxt(f, data, fmt=fmt, header=header)

                    else:
                        header_list = ['']
                        data = np.hstack(
                            (self.energy.reshape(-1, 1), dos[elem][ang]))
                        with open(elem_file_dir/f"{keyname}-{elem}_{ang}.dat", 'w') as f:
                            header_list.append(
                                f"Partial DOS for {keyname}: {elem}")
                            header_list.append('')
                            for orb in self.orbitals:
                                if orb[keyname] == elem and orb["l"] == l_index:
                                    header_list.append(
                                        f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                            header_list.append('')
                            header_list.append(
                                '\tEnergy'+10*' '+'spin 1'+8*' '+'spin 2')
                            header_list.append('')
                            header = '\n'.join(header_list)
                            np.savetxt(f, data, fmt=fmt, header=header)

    def write(self,
              index: Union[Sequence[int], Dict[int, List[int]],
                           Dict[int, Dict[int, List[int]]]] = [],
              atom_index: Union[Sequence[int], Dict[int, List[int]],
                                Dict[int, Dict[int, List[int]]]] = [],
              species: Union[Sequence[str], Dict[str, List[int]],
                             Dict[str, Dict[int, List[int]]]] = [],
              outdir: PathLike = './'
              ):
        """Write parsed partial dos data to files

        Args:
            index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PDOS of each atom. Defaults to [].
            atom_index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PDOS of each atom with same atom_index. Defaults to [].
            species (Union[Sequence[str], Dict[str, List[int]], Dict[str, Dict[int, List[int]]]], optional): extract PDOS of each atom with same species. Defaults to [].
            outdir (PathLike, optional): directory of parsed PDOS files. Defaults to './'.
        """
        if index:
            self._write(index, keyname='index', outdir=outdir)
        if atom_index:
            self._write(atom_index, keyname='atom_index', outdir=outdir)
        if species:
            self._write(species, keyname='species', outdir=outdir)

    def _shift_energy(self, efermi: float = 0, shift: bool = False, prec: float = 0.01):
        tdos = self._all_sum()
        if shift:
            vb, cb = self.set_vcband(
                energy_minus_efermi(self.energy, efermi), tdos, prec)
            self.info(vb, cb)
            energy_f = np.concatenate((vb.band, cb.band))
        else:
            energy_f = energy_minus_efermi(self.energy, efermi)

        return energy_f, tdos

    def _parial_plot(self,
                     fig: Figure,
                     ax: Union[axes.Axes, Sequence[axes.Axes]],
                     species: Union[Sequence[Any], Dict[Any, List[int]],
                                    Dict[Any, Dict[int, List[int]]]] = [],
                     efermi: float = 0,
                     energy_range: Sequence[float] = [],
                     dos_range: Sequence[float] = [],
                     shift: bool = False,
                     prec: float = 0.01,
                     keyname: str = '',
                     **kwargs):
        """Plot parsed partial dos data

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            dos_range (Sequence[float], optional): dos range for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            prec (float, optional): precision for treating dos as zero. Defaults to 0.01.
            keyname (str, optional): the keyword that extracts the PDOS. Defaults to ''.

        Returns:
            DOSPlot object: for manually plotting picture with dosplot.ax 
        """
        if isinstance(ax, axes.Axes):
            ax = [ax]

        dos, totnum = parse_projected_data(self.orbitals, species, keyname)
        energy_f, tdos = self._shift_energy(efermi, shift, prec)

        if not species:
            dosplot = DOSPlot(fig, ax, **kwargs)
            dosplot.ax = self._plot(dosplot, energy_f, tdos, "TDOS")
            if "notes" in dosplot.plot_params.keys():
                dosplot._set_figure(energy_range, dos_range,
                                    notes=dosplot.plot_params["notes"])
            else:
                dosplot._set_figure(energy_range, dos_range)

            return dosplot

        if isinstance(species, (list, tuple)):
            dosplot = DOSPlot(fig, ax, **kwargs)
            if "xlabel_params" in dosplot.plot_params.keys():
                dosplot.ax.set_xlabel("Energy(eV)", **
                                      dosplot.plot_params["xlabel_params"])
            else:
                dosplot.ax.set_xlabel("Energy(eV)", size=25)
            dosplot.ax = self._plot(dosplot, energy_f, tdos, "TDOS")
            for elem in dos.keys():
                dosplot.ax = self._plot(dosplot, energy_f, dos[elem], elem)
            if "notes" in dosplot.plot_params.keys():
                dosplot._set_figure(energy_range, dos_range,
                                    notes=dosplot.plot_params["notes"])
            else:
                dosplot._set_figure(energy_range, dos_range)

            return dosplot

        elif isinstance(species, dict):
            dosplots = []
            assert len(ax) >= len(
                dos.keys()), "There must be enough `axes` to plot."
            for i, elem in enumerate(dos.keys()):
                dosplot = DOSPlot(fig, ax[i], **kwargs)
                for ang in dos[elem].keys():
                    l_index = int(ang)
                    if isinstance(dos[elem][ang], dict):
                        for mag in dos[elem][ang].keys():
                            m_index = int(mag)
                            dosplot.ax = self._plot(
                                dosplot, energy_f, dos[elem][ang][mag], f"{elem}-{get_angular_momentum_name(l_index, m_index)}")

                    else:
                        dosplot.ax = self._plot(
                            dosplot, energy_f, dos[elem][ang], f"{elem}-{get_angular_momentum_label(l_index)}")
                if "notes" in dosplot.plot_params.keys():
                    dosplot._set_figure(energy_range, dos_range,
                                        notes=dosplot.plot_params["notes"][i])
                else:
                    dosplot._set_figure(energy_range, dos_range)
                dosplots.append(dosplot)
            if "xlabel_params" in dosplot.plot_params.keys():
                dosplot.ax.set_xlabel("Energy(eV)", **
                                      dosplot.plot_params["xlabel_params"])
            else:
                dosplot.ax.set_xlabel("Energy(eV)", size=25)

            return dosplots

    def plot(self,
             fig: Figure,
             ax: Union[axes.Axes, Sequence[axes.Axes]],
             index: Union[Sequence[int], Dict[int, List[int]],
                          Dict[int, Dict[int, List[int]]]] = [],
             atom_index: Union[Sequence[int], Dict[int, List[int]],
                               Dict[int, Dict[int, List[int]]]] = [],
             species: Union[Sequence[str], Dict[str, List[int]],
                            Dict[str, Dict[int, List[int]]]] = [],
             efermi: float = 0,
             energy_range: Sequence[float] = [],
             dos_range: Sequence[float] = [],
             shift: bool = False,
             prec: float = 0.01,
             **kwargs):
        """Plot parsed partial dos data

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PDOS of each atom. Defaults to [].
            atom_index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PDOS of each atom with same atom_index. Defaults to [].
            species (Union[Sequence[str], Dict[str, List[int]], Dict[str, Dict[int, List[int]]]], optional): extract PDOS of each atom with same species. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            dos_range (Sequence[float], optional): dos range for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            prec (float, optional): precision for treating dos as zero. Defaults to 0.01.

        Returns:
            DOSPlot object: for manually plotting picture with dosplot.ax 
        """

        if not index and not atom_index and not species:
            dosplot = self._parial_plot(fig=fig, ax=ax, species=[
            ], efermi=efermi, energy_range=energy_range, dos_range=dos_range, shift=shift, prec=prec, keyname='', **kwargs)
        if index:
            dosplot = self._parial_plot(fig=fig, ax=ax, species=index, efermi=efermi, energy_range=energy_range,
                                        dos_range=dos_range, shift=shift, prec=prec, keyname='index', **kwargs)
        if atom_index:
            dosplot = self._parial_plot(fig=fig, ax=ax, species=atom_index, efermi=efermi, energy_range=energy_range,
                                        dos_range=dos_range, shift=shift, prec=prec, keyname='atom_index', **kwargs)
        if species:
            dosplot = self._parial_plot(fig=fig, ax=ax, species=species, efermi=efermi, energy_range=energy_range,
                                        dos_range=dos_range, shift=shift, prec=prec, keyname='species', **kwargs)

        return dosplot


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #tdosfile = r"C:\Users\YY.Ji\Desktop\TDOS"
    #tdos = TDOS(tdosfile)
    #fig, ax = plt.subplots()
    #energy_range = [-6, 10]
    #dos_range = [0, 5]
    # dosplots = tdos.plot(fig, ax, efermi=5, shift=True,
    #                     energy_range=energy_range, dos_range=dos_range, notes={'s': '(a)'})
    # fig.savefig("tdos.png")

    pdosfile = r"../examples/Si/PDOS"
    pdos = PDOS(pdosfile)
    #species = {"Ag": [2], "Cl": [1], "In": [0]}
    atom_index = {1: {1: [0, 1]}}
    fig, ax = plt.subplots(1, 1, sharex=True)
    energy_range = [-5, 7]
    efermi = 6.585653952007503
    dos_range = [0, 5]

    # if you want to specify `species` or `index`, you need to
    # set `species=species` or `index=index` in the following two functions
    dosplots = pdos.plot(fig, ax, atom_index=atom_index, efermi=efermi, shift=True,
                         energy_range=energy_range, dos_range=dos_range, notes=[{'s': '(a)'}])
    fig.savefig("pdos.png")
    pdos.write(atom_index=atom_index)
