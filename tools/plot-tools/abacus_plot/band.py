'''
Date: 2021-12-29 10:27:01
LastEditors: jiyuyang
LastEditTime: 2022-08-02 11:36:22
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
import numpy as np
from os import PathLike
from typing import Sequence, Tuple, Union, Any, List, Dict
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap
from matplotlib.collections import LineCollection
from matplotlib import axes
import matplotlib.pyplot as plt
from pathlib import Path

from abacus_plot.utils import energy_minus_efermi, list_elem2str, read_kpt, remove_empty, parse_projected_data, handle_data, get_angular_momentum_name, get_angular_momentum_label


class Band:
    """Parse Bands data"""

    def __init__(self, bandfile: Union[PathLike, Sequence[PathLike]] = None, kptfile: PathLike = '', old_ver=False) -> None:
        self.bandfile = bandfile
        self.old_ver = old_ver
        if isinstance(bandfile, list) or isinstance(bandfile, tuple):
            self.energy = []
            for file in self.bandfile:
                self.k_index, e, self.k_lines = self.read(file, self.old_ver)
                self.energy.append(e)
        else:
            self.k_index, self.energy, self.k_lines = self.read(
                self.bandfile, self.old_ver)
        self.kptfile = kptfile
        self.kpt = None
        if self.kptfile:
            self.kpt = read_kpt(kptfile)
        self.k_index = list(map(int, self.k_index))
        if self.kpt:
            self._kzip = self.kpt.label_special_k
        else:
            self._kzip = self.k_index

    @classmethod
    def read(cls, filename: PathLike, old_ver=False):
        """Read band data file and return k-points and energy

        :params filename: string of band data file
        """

        z = None
        if old_ver:
            data = np.loadtxt(filename, dtype=float)
            X, y = np.split(data, (1, ), axis=1)
            x = X.flatten()
        else:
            data = np.loadtxt(filename, dtype=float)
            #X, Z, y, _ = np.split(data, [1, 2, data.shape[-1]], axis=1)
            #x = X.flatten()
            #z = Z.flatten()
            x = data[:, 0]     # k-index
            z = data[:, 1]     # k-points
            y = data[:, 2:]    # band

        return x, y, z

    @classmethod
    def direct_bandgap(cls, vb: namedtuple, cb: namedtuple, klength: int):
        """Calculate direct band gap"""

        gap_list = []
        i_index = []
        for i in range(klength):
            gap_list.append(np.min(cb.band[:, i])-np.max(vb.band[:, i]))
            i_index.append(i)
        dgap = np.min(gap_list)

        return dgap, i_index[np.argmin(gap_list)]

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap

    @classmethod
    def band_type(cls, vb: namedtuple, cb: namedtuple):
        vbm_x, cbm_x = vb.k_index, cb.k_index
        longone, shortone = (vbm_x, cbm_x) if len(
            vbm_x) >= len(cbm_x) else (cbm_x, vbm_x)
        for i in shortone:
            if i in longone:
                btype = "Direct"
            else:
                btype = "Indirect"

        return btype

    @classmethod
    def info(cls, kpath: Sequence, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path
        :params energy: band energy after subtracting the Fermi level
        """

        gap = cls.bandgap(vb, cb)
        dgap, d_i = cls.direct_bandgap(vb, cb, len(kpath))
        btype = cls.band_type(vb, cb)
        print(
            "--------------------------Band Structure--------------------------", flush=True)
        print(
            f"{'Band character:'.ljust(30)}{btype}", flush=True)
        if btype == "Indirect":
            print(f"{'Direct Band gap(eV):'.ljust(30)}{dgap: .4f}", flush=True)
            print(f"{'Indirect Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)
        elif btype == "Direct":
            print(f"{'Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)
        print(f"{'Band index:'.ljust(30)}{'HOMO'.ljust(10)}{'LUMO'}", flush=True)
        print(
            f"{''.ljust(30)}{str(vb.band_index[-1]).ljust(10)}{str(cb.band_index[0])}", flush=True)
        print(f"{'Eigenvalue of VBM(eV):'.ljust(30)}{vb.value: .4f}", flush=True)
        print(f"{'Eigenvalue of CBM(eV):'.ljust(30)}{cb.value: .4f}", flush=True)
        vbm_k = np.unique(kpath[vb.k_index], axis=0)
        cbm_k = np.unique(kpath[cb.k_index], axis=0)
        print(
            f"{'Location of VBM'.ljust(30)}{' '.join(list_elem2str(vbm_k[0]))}", flush=True)
        for i, j in enumerate(vbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)
        print(
            f"{'Location of CBM'.ljust(30)}{' '.join(list_elem2str(cbm_k[0]))}", flush=True)
        for i, j in enumerate(cbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)

    @classmethod
    def set_vcband(cls, energy: Sequence) -> Tuple[namedtuple, namedtuple]:
        """Separate valence and conduct band

        :params energy: band energy after subtracting the Fermi level
        """

        e_T = energy.T
        num_gt_Ef = (e_T > 0).sum(axis=1)

        Band = namedtuple('Band', ['band_index', 'band', 'value', 'k_index'])

        # valance band
        band_vbm_index = np.where(num_gt_Ef == 0)[0]
        band_vbm = e_T[band_vbm_index]
        evbm = np.max(band_vbm)
        k_vbm_index = np.where(band_vbm == evbm)[1]
        vb = Band(band_vbm_index, band_vbm-evbm, evbm-evbm, k_vbm_index)

        # conduct band
        band_cbm_index = np.where(num_gt_Ef != 0)[0]
        band_cbm = e_T[band_cbm_index]
        ecbm = np.min(band_cbm)
        k_cbm_index = np.where(band_cbm == ecbm)[1]
        cb = Band(band_cbm_index, band_cbm-evbm, ecbm-evbm, k_cbm_index)

        return vb, cb

    def _shift_energy(self, energy, efermi: float = 0, shift: bool = False):
        energy = energy_minus_efermi(energy, efermi)
        if shift:
            vb, cb = self.set_vcband(energy)
            refine_E = np.vstack((vb.band, cb.band)).T
            self.info(self.kpt.full_kpath, vb, cb)
        else:
            refine_E = energy

        return refine_E

    @classmethod
    def plot_data(cls,
                  fig: Figure,
                  ax: axes.Axes,
                  x: Sequence,
                  y: Sequence,
                  kzip: Sequence,
                  efermi: float = 0,
                  energy_range: Sequence[float] = [],
                  **kwargs):
        """Plot band structure

        :params x, y: x-axis and y-axis coordinates
        :params kzip: special k-points label and its k-index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        """

        bandplot = BandPlot(fig, ax, **kwargs)
        if not bandplot._color:
            bandplot._color = 'black'

        kpoints, energy = x, y
        energy = energy_minus_efermi(energy, efermi)

        bandplot.ax.plot(kpoints, energy, lw=bandplot._lw, color=bandplot._color,
                         label=bandplot._label, linestyle=bandplot._linestyle)
        bandplot._set_figure(kzip, energy_range)

    def plot(self,
             fig: Figure,
             ax: axes.Axes,
             efermi: Union[float, Sequence[float]] = [],
             energy_range: Sequence[float] = [],
             shift: bool = True,
             **kwargs):
        """Plot more than two band structures using data file

        :params efermi: Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        """

        bandplot = BandPlot(fig, ax, **kwargs)

        if isinstance(self.bandfile, list):
            nums = len(self.bandfile)
            if not efermi:
                efermi = [0.0 for i in range(nums)]
            if not kwargs.pop('label', None):
                bandplot._label = ['' for i in range(nums)]
            if not kwargs.pop('color', None):
                bandplot._color = ['black' for i in range(nums)]
            if not kwargs.pop('linestyle', None):
                bandplot._linestyle = ['solid' for i in range(nums)]

            for i, band in enumerate(self.energy):
                band = self._shift_energy(band, efermi[i], shift)
                bandplot.ax.plot(self.k_index, band,
                                 lw=bandplot._lw, color=bandplot._color[i], label=bandplot._label[i], linestyle=bandplot._linestyle[i])

        else:
            if not efermi:
                efermi = 0.0

            band = self._shift_energy(self.energy, efermi, shift)
            bandplot.ax.plot(self.k_index, band,
                             lw=bandplot._lw, color=bandplot._color, label=bandplot._label, linestyle=bandplot._linestyle)

        bandplot._set_figure(self._kzip, energy_range)

        return bandplot


class BandPlot:
    """Plot band structure"""

    def __init__(self, fig: Figure, ax: axes.Axes, point_to_line=False, **kwargs) -> None:
        self.fig = fig
        self.ax = ax
        self._point_to_line = point_to_line
        self._lw = kwargs.pop('lw', 2)
        self._bwidth = kwargs.pop('bwdith', 3)
        self._label = kwargs.pop('label', None)
        self._color = kwargs.pop('color', 'black')
        self._linestyle = kwargs.pop('linestyle', 'solid')
        self.plot_params = kwargs

    def _set_figure(self, kzip, range: Sequence):
        """set figure and axes for plotting

        :params kzip: dict of label of points of x-axis and its index in data file. Range of x-axis based on kzip.value()
        :params range: range of y-axis
        """

        keys = []
        values = []
        for t in kzip:
            if isinstance(t, tuple) or isinstance(t, list):
                keys.append(t[0])
                values.append(t[1])
            else:
                keys.append('')
                values.append(int(t))

        # x-axis
        self.ax.set_xticks(values)
        self.ax.set_xticklabels(keys)
        self.ax.set_xlim(values[0], values[-1])
        if "xlabel_params" in self.plot_params.keys():
            self.ax.set_xlabel(
                "Wave Vector", **self.plot_params["xlabel_params"])
        else:
            self.ax.set_xlabel("Wave Vector", size=25)

        # y-axis
        if range:
            self.ax.set_ylim(range[0], range[1])
        if "ylabel_params" in self.plot_params.keys():
            self.ax.set_ylabel(
                "Energy(eV)", **self.plot_params["ylabel_params"])
        else:
            self.ax.set_ylabel("Energy(eV)", size=25)

        # notes
        if "notes" in self.plot_params.keys():
            from matplotlib.offsetbox import AnchoredText
            if "s" in self.plot_params["notes"].keys() and len(self.plot_params["notes"].keys()) == 1:
                self.ax.add_artist(AnchoredText(self.plot_params["notes"]["s"], loc='upper left', prop=dict(size=25),
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

        # guides
        if '' not in keys:
            if "grid_params" in self.plot_params.keys():
                self.ax.grid(axis='x', **self.plot_params["grid_params"])
            else:
                self.ax.grid(axis='x', lw=1.2)
        if "hline_params" in self.plot_params.keys():
            self.ax.axhline(0, **self.plot_params["hline_params"])
        else:
            self.ax.axhline(0, linestyle="--", c='b', lw=1.0)

        if self._label:
            handles, labels = self.ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if "legend_prop" in self.plot_params.keys():
                self.ax.legend(by_label.values(), by_label.keys(),
                               prop=self.plot_params["legend_prop"])
            else:
                self.ax.legend(by_label.values(),
                               by_label.keys(), prop={'size': 15})

        if self._point_to_line:
            handles = []
            for c, l in zip(self._color, self._label):
                handles.append(Line2D([0], [0], color=c, label=l))
            if "legend_prop" in self.plot_params.keys():
                self.ax.legend(handles=handles,
                               prop=self.plot_params["legend_prop"])
            else:
                self.ax.legend(handles=handles, prop={'size': 15})

    def _color_to_alpha_cmap(self, color):
        cmap = LinearSegmentedColormap.from_list("chaos", ["white", color])
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:, -1] = np.linspace(0, 1, cmap.N)  # this adds alpha
        my_cmap = ListedColormap(my_cmap)
        return my_cmap


class PBand(Band):
    def __init__(self, bandfile: Union[PathLike, Sequence[PathLike]] = None, kptfile: str = '') -> None:
        self.bandfile = bandfile
        if isinstance(bandfile, list) or isinstance(bandfile, tuple):
            self.energy = []
            self.orbitals = []
            for file in self.bandfile:
                self.nspin, self.norbitals, self.eunit, self.nbands, self.nkpoints, self.k_index, e, orb = self.read(
                    file)
                self._check_energy(e)
                self.energy.append(e)
                self.orbitals.append(orb)
        else:
            self.nspin, self.norbitals, self.eunit, self.nbands, self.nkpoints, self.k_index, self.energy, self.orbitals = self.read(
                self.bandfile)
            self._check_energy(self.energy)
        self.energy = np.asarray(self.energy)
        self.kptfile = kptfile
        self.kpt = None
        if self.kptfile:
            self.kpt = read_kpt(kptfile)
        self.k_index = list(map(int, self.k_index))
        if self.kpt:
            self._kzip = self.kpt.label_special_k
        else:
            self._kzip = self.k_index

        self._check_weights(self.weights)

    def _check_energy(self, energy):
        assert energy.shape[0] == self.nkpoints, "The dimension of band structure dismatches with the number of k-points."
        assert energy.shape[1] == self.nbands, "The dimension of band structure dismatches with the number of bands."

    def _check_weights(self, weights: np.ndarray, prec=1e-5):
        assert weights.shape[0] == self.norbitals, "The dimension of weights dismatches with the number of orbitals."
        assert weights.shape[1] == self.nkpoints, "The dimension of weights dismatches with the number of k-points."
        assert weights.shape[2] == self.nbands, "The dimension of weights dismatches with the number of bands."
        one_mat = np.ones((self.nkpoints, self.nbands))
        assert (np.abs(weights.sum(axis=0)-one_mat) < prec).all(
        ), f"np.abs(weights.sum(axis=0)-np.ones(({self.nkpoints}, {self.nbands}))) < {prec}"

    @property
    def weights(self):
        data = np.empty((self.norbitals, self.nkpoints, self.nbands))
        for i, orb in enumerate(self.orbitals):
            data[i] = orb['data']
        return data

    @classmethod
    def read(cls, filename: PathLike):
        """Read projected band data file and return k-points, energy and Mulliken weights

        :params bandfile: string of projected band data file
        """

        from lxml import etree
        pbanddata = etree.parse(filename)
        root = pbanddata.getroot()
        nspin = int(root.xpath('//nspin')[0].text.replace(' ', ''))
        norbitals = int(root.xpath('//norbitals')
                        [0].text.replace(' ', ''))
        eunit = root.xpath('//band_structure/@units')[0].replace(' ', '')
        nbands = int(root.xpath('//band_structure/@nbands')
                     [0].replace(' ', ''))
        nkpoints = int(root.xpath('//band_structure/@nkpoints')
                       [0].replace(' ', ''))
        k_index = np.arange(nkpoints)
        energy = root.xpath('//band_structure')[0].text.split('\n')
        energy = handle_data(energy)
        remove_empty(energy)
        energy = np.asarray(energy, dtype=float)

        orbitals = []
        for i in range(norbitals):
            orb = OrderedDict()
            o_index_str = root.xpath(
                '//orbital/@index')[i]
            orb['index'] = int(o_index_str.replace(' ', ''))
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
            orbitals.append(orb)

        return nspin, norbitals, eunit, nbands, nkpoints, k_index, energy, orbitals

    def _write(self, species: Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], keyname='', file_dir: PathLike = ''):
        """Write parsed projected bands data to files

        Args:
            orbital (dict): parsed data
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            keyname (str): the keyword that extracts the PBAND. Allowed values: 'index', 'atom_index', 'species'
        """

        band, totnum = parse_projected_data(self.orbitals, species, keyname)

        if isinstance(species, (list, tuple)):
            for elem in band.keys():
                header_list = ['']
                with open(file_dir/f"{keyname}-{elem}.dat", 'w') as f:
                    header_list.append(
                        f"Projected band structure for {keyname}: {elem}")
                    header_list.append('')
                    header_list.append(
                        f'\tNumber of k-points: {self.nkpoints}')
                    header_list.append(f'\tNumber of bands: {self.nbands}')
                    header_list.append('')
                    for orb in self.orbitals:
                        if orb[keyname] == elem:
                            header_list.append(
                                f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                    header_list.append('')
                    header_list.append(
                        f'Data shape: ({self.nkpoints}, {self.nbands})')
                    header_list.append('')
                    header = '\n'.join(header_list)
                    np.savetxt(f, band[elem], header=header)

        elif isinstance(species, dict):
            for elem in band.keys():
                elem_file_dir = file_dir/f"{keyname}-{elem}"
                elem_file_dir.mkdir(exist_ok=True)
                for ang in band[elem].keys():
                    l_index = int(ang)
                    if isinstance(band[elem][ang], dict):
                        for mag in band[elem][ang].keys():
                            header_list = ['']
                            m_index = int(mag)
                            with open(elem_file_dir/f"{keyname}-{elem}_{ang}_{mag}.dat", 'w') as f:
                                header_list.append(
                                    f"Projected band structure for {keyname}: {elem}")
                                header_list.append('')
                                header_list.append(
                                    f'\tNumber of k-points: {self.nkpoints}')
                                header_list.append(
                                    f'\tNumber of bands: {self.nbands}')
                                header_list.append('')
                                for orb in self.orbitals:
                                    if orb[keyname] == elem and orb["l"] == l_index and orb["m"] == m_index:
                                        header_list.append(
                                            f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                                header_list.append('')
                                header_list.append(
                                    f'Data shape: ({self.nkpoints}, {self.nbands})')
                                header_list.append('')
                                header = '\n'.join(header_list)
                                np.savetxt(f, band[elem][ang]
                                           [mag], header=header)

                    else:
                        header_list = ['']
                        with open(elem_file_dir/f"{keyname}-{elem}_{ang}.dat", 'w') as f:
                            header_list.append(
                                f"Projected band structure for {keyname}: {elem}")
                            header_list.append('')
                            header_list.append(
                                f'\tNumber of k-points: {self.nkpoints}')
                            header_list.append(
                                f'\tNumber of bands: {self.nbands}')
                            header_list.append('')
                            for orb in self.orbitals:
                                if orb[keyname] == elem and orb["l"] == l_index:
                                    header_list.append(
                                        f"\tAdd data for index ={orb['index']:4d}, atom_index ={orb['atom_index']:4d}, element ={orb['species']:4s},  l,m,z={orb['l']:3d}, {orb['m']:3d}, {orb['z']:3d}")
                            header_list.append('')
                            header_list.append(
                                f'Data shape: ({self.nkpoints}, {self.nbands})')
                            header_list.append('')
                            header = '\n'.join(header_list)
                            np.savetxt(f, band[elem][ang], header=header)

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

        if isinstance(self.bandfile, list):
            for i in range(len(self.bandfile)):
                file_dir = Path(f"{outdir}", f"PBAND{i}_FILE")
                file_dir.mkdir(exist_ok=True)
                if index:
                    self._write(index, keyname='index',
                                file_dir=file_dir)
                if atom_index:
                    self._write(atom_index, keyname='atom_index',
                                file_dir=file_dir)
                if species:
                    self._write(species, keyname='species',
                                file_dir=file_dir)

        else:
            file_dir = Path(f"{outdir}", f"PBAND{1}_FILE")
            file_dir.mkdir(exist_ok=True)
            if index:
                self._write(index, keyname='index',
                            file_dir=file_dir)
            if atom_index:
                self._write(atom_index, keyname='atom_index',
                            file_dir=file_dir)
            if species:
                self._write(species, keyname='species',
                            file_dir=file_dir)

    def _plot(self,
              fig: Figure,
              ax: axes.Axes,
              energy: np.ndarray,
              species: Union[Sequence[Any], Dict[Any, List[int]],
                             Dict[Any, Dict[int, List[int]]]] = [],
              efermi: float = 0,
              energy_range: Sequence[float] = [],
              shift: bool = False,
              keyname: str = '',
              outdir: PathLike = './',
              out_index: int = 1,
              cmap='jet',
              **kwargs):
        """Plot parsed projected bands data

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            keyname (str, optional): the keyword that extracts the PBANDS. Defaults to ''.

        Returns:
            BandPlot object: for manually plotting picture with bandplot.ax 
        """

        def _seg_plot(bandplot, lc, file_dir, name):
            cbar = bandplot.fig.colorbar(lc, ax=bandplot.ax)
            bandplot._set_figure(self._kzip, energy_range)
            bandplot.fig.savefig(
                file_dir/f'{keyname}-{bandplot._label}.pdf', dpi=400)
            cbar.remove()
            plt.cla()

            return bandplot

        wei, totnum = parse_projected_data(self.orbitals, species, keyname)
        energy = self._shift_energy(energy, efermi, shift)
        file_dir = Path(f"{outdir}", f"PBAND{out_index}_FIG")
        file_dir.mkdir(exist_ok=True)

        if not species:
            bandplot = BandPlot(fig, ax, **kwargs)
            bandplot = super().plot(fig, ax, efermi, energy_range, shift, **kwargs)
            bandplot._set_figure(self._kzip, energy_range)

            return bandplot

        if isinstance(species, (list, tuple)):
            bandplots = []
            for i, elem in enumerate(wei.keys()):
                bandplot = BandPlot(fig, ax, **kwargs)
                bandplot._label = elem
                for ib in range(self.nbands):
                    points = np.array(
                        (self.k_index, energy[0:, ib])).T.reshape(-1, 1, 2)
                    segments = np.concatenate(
                        [points[:-1], points[1:]], axis=1)
                    norm = Normalize(
                        vmin=wei[elem][0:, ib].min(), vmax=wei[elem][0:, ib].max())
                    lc = LineCollection(
                        segments, cmap=plt.get_cmap(cmap), norm=norm)
                    lc.set_array(wei[elem][0:, ib])
                    lc.set_label(bandplot._label)
                    bandplot.ax.add_collection(lc)

                _seg_plot(bandplot, lc, file_dir, name=f'{elem}')
                bandplots.append(bandplot)
            return bandplots

        elif isinstance(species, dict):
            bandplots = []
            for i, elem in enumerate(wei.keys()):
                elem_file_dir = file_dir/f"{keyname}-{elem}"
                elem_file_dir.mkdir(exist_ok=True)
                for ang in wei[elem].keys():
                    l_index = int(ang)
                    if isinstance(wei[elem][ang], dict):
                        for mag in wei[elem][ang].keys():
                            bandplot = BandPlot(fig, ax, **kwargs)
                            m_index = int(mag)
                            bandplot._label = f"{elem}-{get_angular_momentum_name(l_index, m_index)}"
                            for ib in range(self.nbands):
                                points = np.array(
                                    (self.k_index, energy[0:, ib])).T.reshape(-1, 1, 2)
                                segments = np.concatenate(
                                    [points[:-1], points[1:]], axis=1)
                                norm = Normalize(vmin=wei[elem][ang][mag][0:, ib].min(
                                ), vmax=wei[elem][ang][mag][0:, ib].max())
                                lc = LineCollection(
                                    segments, cmap=plt.get_cmap(cmap), norm=norm)
                                lc.set_array(wei[elem][ang][mag][0:, ib])
                                lc.set_label(bandplot._label)
                                bandplot.ax.add_collection(lc)

                            _seg_plot(bandplot, lc, elem_file_dir,
                                      name=f'{elem}_{ang}_{mag}')
                            bandplots.append(bandplot)

                    else:
                        bandplot = BandPlot(fig, ax, **kwargs)
                        bandplot._label = f"{elem}-{get_angular_momentum_label(l_index)}"
                        for ib in range(self.nbands):
                            points = np.array(
                                (self.k_index, energy[0:, ib])).T.reshape(-1, 1, 2)
                            segments = np.concatenate(
                                [points[:-1], points[1:]], axis=1)
                            norm = Normalize(vmin=wei[elem][ang][0:, ib].min(
                            ), vmax=wei[elem][ang][0:, ib].max())
                            lc = LineCollection(
                                segments, cmap=plt.get_cmap(cmap), norm=norm)
                            lc.set_array(wei[elem][ang][0:, ib])
                            lc.set_label(bandplot._label)
                            bandplot.ax.add_collection(lc)

                        _seg_plot(bandplot, lc, elem_file_dir,
                                  name=f'{elem}_{ang}')
                        bandplots.append(bandplot)

            return bandplots

        plt.clf()

    def plot(self,
             fig: Figure,
             ax: Union[axes.Axes, Sequence[axes.Axes]],
             index: Union[Sequence[int], Dict[int, List[int]],
                          Dict[int, Dict[int, List[int]]]] = [],
             atom_index: Union[Sequence[int], Dict[int, List[int]],
                               Dict[int, Dict[int, List[int]]]] = [],
             species: Union[Sequence[str], Dict[str, List[int]],
                            Dict[str, Dict[int, List[int]]]] = [],
             efermi: Union[float, Sequence[float]] = [],
             energy_range: Sequence[float] = [],
             shift: bool = False,
             outdir: PathLike = './',
             cmapname='jet',
             **kwargs):
        """Plot parsed projected band data

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PBAND of each atom. Defaults to [].
            atom_index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PBAND of each atom with same atom_index. Defaults to [].
            species (Union[Sequence[str], Dict[str, List[int]], Dict[str, Dict[int, List[int]]]], optional): extract PBAND of each atom with same species. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            outdir (PathLike): Default: './'
            cmapname (str): Default: 'jet'

        Returns:
            BandPlot object: for manually plotting picture with bandplot.ax 
        """

        if isinstance(self.bandfile, list):
            nums = len(self.bandfile)
            if not efermi:
                efermi = [0.0 for i in range(nums)]
            _linestyle = kwargs.pop(
                'linestyle', ['solid' for i in range(nums)])

            for i, band in enumerate(self.energy):
                if not index and not atom_index and not species:
                    bandplot = self._plot(fig=fig, ax=ax, energy=band, species=[
                    ], efermi=efermi[i], energy_range=energy_range, shift=shift, keyname='', linestyle=_linestyle[i], outdir=outdir, out_index=i, cmapname=cmapname, **kwargs)
                if index:
                    bandplot = self._plot(fig=fig, ax=ax, energy=band, species=index, efermi=efermi[i],
                                          energy_range=energy_range, shift=shift, keyname='index', linestyle=_linestyle[i], outdir=outdir, out_index=i, cmapname=cmapname, **kwargs)
                if atom_index:
                    bandplot = self._plot(fig=fig, ax=ax, energy=band, species=atom_index, efermi=efermi[i],
                                          energy_range=energy_range, shift=shift, keyname='atom_index', linestyle=_linestyle[i], outdir=outdir, out_index=i, cmapname=cmapname, **kwargs)
                if species:
                    bandplot = self._plot(fig=fig, ax=ax, energy=band, species=species, efermi=efermi[i],
                                          energy_range=energy_range, shift=shift, keyname='species', linestyle=_linestyle[i], outdir=outdir, out_index=i, cmapname=cmapname, **kwargs)

        else:
            if not index and not atom_index and not species:
                bandplot = self._plot(fig=fig, ax=ax, energy=self.energy, species=[
                ], efermi=efermi, energy_range=energy_range, shift=shift, keyname='', outdir=outdir, out_index=1, **kwargs)
            if index:
                bandplot = self._plot(fig=fig, ax=ax, energy=self.energy, species=index, efermi=efermi,
                                      energy_range=energy_range, shift=shift, keyname='index', outdir=outdir, out_index=1, **kwargs)
            if atom_index:
                bandplot = self._plot(fig=fig, ax=ax, energy=self.energy, species=atom_index, efermi=efermi,
                                      energy_range=energy_range, shift=shift, keyname='atom_index', outdir=outdir, out_index=1, **kwargs)
            if species:
                bandplot = self._plot(fig=fig, ax=ax, energy=self.energy, species=species, efermi=efermi,
                                      energy_range=energy_range, shift=shift, keyname='species', outdir=outdir, out_index=1, **kwargs)

        return bandplot

    def _plot_contributions(self,
                            fig: Figure,
                            ax: axes.Axes,
                            energy: np.ndarray,
                            species: Union[Sequence[Any], Dict[Any, List[int]],
                                           Dict[Any, Dict[int, List[int]]]] = [],
                            efermi: float = 0,
                            energy_range: Sequence[float] = [],
                            shift: bool = False,
                            keyname: str = '',
                            colors: list = [],
                            scale_width_factor: int = 5,
                            **kwargs):
        """Plot parsed projected bands data of different contributions

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            species (Union[Sequence[Any], Dict[Any, List[int]], Dict[Any, Dict[int, List[int]]]], optional): list of atomic species(index or atom index) or dict of atomic species(index or atom index) and its angular momentum list. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            keyname (str, optional): the keyword that extracts the PBANDS. Defaults to ''.

        Returns:
            BandPlot object: for manually plotting picture with bandplot.ax 
        """

        wei, totnum = parse_projected_data(self.orbitals, species, keyname)
        energy = self._shift_energy(energy, efermi, shift)

        if not species:
            bandplot = BandPlot(fig, ax, **kwargs)
            bandplot = super().plot(fig, ax, efermi, energy_range, shift, **kwargs)
            bandplot._set_figure(self._kzip, energy_range)

            return bandplot

        whole_data_parsed = []
        whole_label_parsed = []
        if isinstance(species, (list, tuple)):
            for i, elem in enumerate(wei.keys()):
                whole_label_parsed.append(elem)
                whole_data_parsed.append(wei[elem])

        elif isinstance(species, dict):
            for i, elem in enumerate(wei.keys()):
                for ang in wei[elem].keys():
                    l_index = int(ang)
                    if isinstance(wei[elem][ang], dict):
                        for mag in wei[elem][ang].keys():
                            m_index = int(mag)
                            whole_label_parsed.append(
                                f"{elem}-{get_angular_momentum_name(l_index, m_index)}")
                            whole_data_parsed.append(wei[elem][ang][mag])

                    else:
                        whole_label_parsed.append(
                            f"{elem}-{get_angular_momentum_label(l_index)}")
                        whole_data_parsed.append(wei[elem][ang])

        if len(colors) == 0:
            cmap = plt.cm.get_cmap("tab10")
            colors = [cmap(c)
                      for c in np.linspace(0, 1, len(whole_label_parsed))]

        bandplot = BandPlot(
            fig, ax, color=colors, label=whole_label_parsed, point_to_line=True, **kwargs)

        norm = Normalize(vmin=0, vmax=1)
        cmaps = [bandplot._color_to_alpha_cmap(c) for c in colors]
        swidth = np.array(whole_data_parsed)*scale_width_factor
        for ib in range(self.nbands):
            for i, con in enumerate(swidth):
                bandplot.ax.scatter(
                    self.k_index, energy[:, ib], c=con[:, ib], s=con[:, ib], norm=norm, cmap=cmaps[i], label=whole_label_parsed[i])

        # clb = plt.colorbar(
        #     plt.cm.ScalarMappable(norm=norm, cmap=cmaps), ax=bandplot.ax
        # )
        # clb.set_ticks(np.linspace(0, 1, len(whole_label_parsed)))
        # clb.set_ticklabels(whole_label_parsed)
        bandplot._set_figure(self._kzip, energy_range)

        return bandplot

    def plot_contributions(self,
                           fig: Figure,
                           ax: Union[axes.Axes, Sequence[axes.Axes]],
                           index: Union[Sequence[int], Dict[int, List[int]],
                                        Dict[int, Dict[int, List[int]]]] = [],
                           atom_index: Union[Sequence[int], Dict[int, List[int]],
                                             Dict[int, Dict[int, List[int]]]] = [],
                           species: Union[Sequence[str], Dict[str, List[int]],
                                          Dict[str, Dict[int, List[int]]]] = [],
                           efermi: Union[float, Sequence[float]] = [],
                           energy_range: Sequence[float] = [],
                           shift: bool = False,
                           colors: list = [],
                           **kwargs):
        """Plot parsed projected band data of different contributions

        Args:
            fig (Figure): object of matplotlib.figure.Figure
            ax (Union[axes.Axes, Sequence[axes.Axes]]): object of matplotlib.axes.Axes or a list of this objects
            index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PBAND of each atom. Defaults to [].
            atom_index (Union[Sequence[int], Dict[int, List[int]], Dict[int, Dict[int, List[int]]]], optional): extract PBAND of each atom with same atom_index. Defaults to [].
            species (Union[Sequence[str], Dict[str, List[int]], Dict[str, Dict[int, List[int]]]], optional): extract PBAND of each atom with same species. Defaults to [].
            efermi (float, optional): fermi level in unit eV. Defaults to 0.
            energy_range (Sequence[float], optional): energy range in unit eV for plotting. Defaults to [].
            shift (bool, optional): if shift energy by fermi level and set the VBM to zero, or not. Defaults to False.
            colors (list, optional): Default:[]

        Returns:
            BandPlot object: for manually plotting picture with bandplot.ax 
        """

        if isinstance(self.bandfile, list):
            nums = len(self.bandfile)
            if not efermi:
                efermi = [0.0 for i in range(nums)]
            _linestyle = kwargs.pop(
                'linestyle', ['solid' for i in range(nums)])

            for i, band in enumerate(self.energy):
                if not index and not atom_index and not species:
                    bandplot = self._plot_contributions(fig=fig, ax=ax, energy=band, species=[
                    ], efermi=efermi[i], energy_range=energy_range, shift=shift, keyname='', linestyle=_linestyle[i], **kwargs)
                if index:
                    bandplot = self._plot_contributions(fig=fig, ax=ax, energy=band, species=index, efermi=efermi[i],
                                                        energy_range=energy_range, shift=shift, keyname='index', linestyle=_linestyle[i], colors=colors, **kwargs)
                if atom_index:
                    bandplot = self._plot_contributions(fig=fig, ax=ax, energy=band, species=atom_index, efermi=efermi[i],
                                                        energy_range=energy_range, shift=shift, keyname='atom_index', linestyle=_linestyle[i], colors=colors, **kwargs)
                if species:
                    bandplot = self._plot_contributions(fig=fig, ax=ax, energy=band, species=species, efermi=efermi[i],
                                                        energy_range=energy_range, shift=shift, keyname='species', linestyle=_linestyle[i], colors=colors, **kwargs)

        else:
            if not index and not atom_index and not species:
                bandplot = self._plot_contributions(fig=fig, ax=ax, energy=self.energy, species=[
                ], efermi=efermi, energy_range=energy_range, shift=shift, keyname='', colors=colors, **kwargs)
            if index:
                bandplot = self._plot_contributions(fig=fig, ax=ax, energy=self.energy, species=index, efermi=efermi,
                                                    energy_range=energy_range, shift=shift, keyname='index', colors=colors, **kwargs)
            if atom_index:
                bandplot = self._plot_contributions(fig=fig, ax=ax, energy=self.energy, species=atom_index, efermi=efermi,
                                                    energy_range=energy_range, shift=shift, keyname='atom_index', colors=colors, **kwargs)
            if species:
                bandplot = self._plot_contributions(fig=fig, ax=ax, energy=self.energy, species=species, efermi=efermi,
                                                    energy_range=energy_range, shift=shift, keyname='species', colors=colors, **kwargs)

        return bandplot


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pathlib import Path
    parent = Path(r"../examples/Si")
    name = "PBANDS_1"
    path = parent/name
    kptfile = parent/'KLINES'
    fig, ax = plt.subplots(figsize=(12, 6))
    energy_range = [-5, 7]
    efermi = 6.585653952007503
    shift = False
    #species = {"Ag": [2], "Cl": [1], "In": [0]}
    atom_index = {1: {1: [0, 1]}}
    pband = PBand(str(path), kptfile)

    # if you want to specify `species` or `index`, you need to
    # set `species=species` or `index=index` in the following two functions

    # 1. write different contributions to files
    # pband.write(atom_index=atom_index)

    # 2. plot different contributions in single picture
    pband.plot_contributions(fig, ax, atom_index=atom_index, efermi=efermi,
                             energy_range=energy_range, shift=shift)
    plt.show()

    # 3. plot different contributions to different pictures with colobar denoting weightes
    # pband.plot(fig, ax, atom_index=atom_index, efermi=efermi,
    #           energy_range=energy_range, shift=shift)       energy_range=energy_range, shift=shift)
