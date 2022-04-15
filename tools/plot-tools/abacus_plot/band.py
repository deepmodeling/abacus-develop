'''
Date: 2021-12-29 10:27:01
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:06:14
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import OrderedDict, namedtuple
import numpy as np
from os import PathLike
from typing import Sequence, Tuple
from matplotlib.figure import Figure
from matplotlib import axes

from abacus_plot.utils import energy_minus_efermi, list_elem2str, read_kpt


class BandPlot:
    """Plot band structure"""

    def __init__(self, fig: Figure, ax: axes.Axes, **kwargs) -> None:
        self.fig = fig
        self.ax = ax
        self._lw = kwargs.pop('lw', 2)
        self._bwidth = kwargs.pop('bwdith', 3)
        self._label = kwargs.pop('label', None)
        self._color = kwargs.pop('color', None)
        self._linestyle = kwargs.pop('linestyle', 'solid')
        self.plot_params = kwargs

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

    @classmethod
    def read(cls, filename: PathLike) -> Tuple[np.ndarray, np.ndarray]:
        """Read band data file and return k-points and energy

        :params filename: string of band data file
        """

        data = np.loadtxt(filename)
        X, y = np.split(data, (1, ), axis=1)
        x = X.flatten()

        return x, y

    def _set_figure(self, index: dict, range: Sequence):
        """set figure and axes for plotting

        :params index: dict of label of points of x-axis and its index in data file. Range of x-axis based on index.value()
        :params range: range of y-axis
        """

        keys = []
        values = []
        for t in index:
            if isinstance(t, tuple):
                keys.append(t[0])
                values.append(t[1])
            elif isinstance(t, (int, float)):
                keys.append('')
                values.append(t)

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

    def plot(self, x: Sequence, y: Sequence, index: Sequence, efermi: float = 0, energy_range: Sequence[float] = []):
        """Plot band structure

        :params x, y: x-axis and y-axis coordinates
        :params index: special k-points label and its index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        """

        if not self._color:
            self._color = 'black'

        kpoints, energy = x, y
        energy = energy_minus_efermi(energy, efermi)

        self.ax.plot(kpoints, energy, lw=self._lw, color=self._color,
                     label=self._label, linestyle=self._linestyle)
        self._set_figure(index, energy_range)

    def singleplot(self, datafile: PathLike, kptfile: str = '', efermi: float = 0, energy_range: Sequence[float] = [], shift: bool = False):
        """Plot band structure using data file

        :params datafile: string of band date file
        :params kptfile: k-point file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        """

        kpt = read_kpt(kptfile)

        if not self._color:
            self._color = 'black'

        kpoints, energy = self.read(datafile)
        energy = energy_minus_efermi(energy, efermi)
        if shift:
            vb, cb = self.set_vcband(energy)
            self.ax.plot(kpoints, np.vstack((vb.band, cb.band)).T,
                         lw=self._lw, color=self._color, label=self._label, linestyle=self._linestyle)
            self.info(kpt.full_kpath, vb, cb)
        else:
            self.ax.plot(kpoints, energy,
                         lw=self._lw, color=self._color, label=self._label, linestyle=self._linestyle)
        index = kpt.label_special_k
        self._set_figure(index, energy_range)

    def multiplot(self, datafile: Sequence[PathLike], kptfile: str = '', efermi: Sequence[float] = [], energy_range: Sequence[float] = [], shift: bool = True):
        """Plot more than two band structures using data file

        :params datafile: list of path of band date file 
        :params kptfile: k-point file
        :params efermi: list of Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params shift: if sets True, it will calculate band gap. This parameter usually is suitable for semiconductor and insulator. Default: False
        """

        kpt = read_kpt(kptfile)

        if not efermi:
            efermi = [0.0 for i in range(len(datafile))]
        if not self._label:
            self._label = ['' for i in range(len(datafile))]
        if not self._color:
            self._color = ['black' for i in range(len(datafile))]
        if not self._linestyle:
            self._linestyle = ['solid' for i in range(len(datafile))]

        emin = -np.inf
        emax = np.inf
        for i, file in enumerate(datafile):
            kpoints, energy = self.read(file)
            if shift:
                vb, cb = self.set_vcband(
                    energy_minus_efermi(energy, efermi[i]))
                energy_min = np.min(vb.band)
                energy_max = np.max(cb.band)
                if energy_min > emin:
                    emin = energy_min
                if energy_max < emax:
                    emax = energy_max

                self.ax.plot(kpoints, np.vstack((vb.band, cb.band)).T,
                             lw=self._lw, color=self._color[i], label=self._label[i], linestyle=self._linestyle[i])
                self.info(kpt.full_kpath, vb, cb)
            else:
                self.ax.plot(kpoints, energy_minus_efermi(energy, efermi[i]),
                             lw=self._lw, color=self._color[i], label=self._label[i], linestyle=self._linestyle[i])

        index = kpt.label_special_k
        self._set_figure(index, energy_range)

    @classmethod
    def bandgap(cls, vb: namedtuple, cb: namedtuple):
        """Calculate band gap"""

        gap = cb.value-vb.value

        return gap

    @classmethod
    def info(cls, kpath: Sequence, vb: namedtuple, cb: namedtuple):
        """Output the information of band structure

        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        def band_type(vbm_x, cbm_x):
            longone, shortone = (vbm_x, cbm_x) if len(
                vbm_x) >= len(cbm_x) else (cbm_x, vbm_x)
            for i in shortone:
                if i in longone:
                    btype = "Direct"
                else:
                    btype = "Indirect"
            return btype

        gap = cls.bandgap(vb, cb)
        print(
            "--------------------------Band Structure--------------------------", flush=True)
        print(
            f"{'Band character:'.ljust(30)}{band_type(vb.k_index, cb.k_index)}", flush=True)
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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pathlib import Path
    parent = Path(r"D:\ustc\TEST\HOIP\double HOIP\result\bond")
    name = "CsAgBiBr"
    path = parent/name
    notes = {'s': '(b)'}
    datafile = [path/"soc.dat", path/"non-soc.dat"]
    kptfile = path/"KPT"
    fig, ax = plt.subplots(figsize=(12, 12))
    label = ["with SOC", "without SOC"]
    color = ["r", "g"]
    linestyle = ["solid", "dashed"]
    band = BandPlot(fig, ax, notes=notes, label=label,
                    color=color, linestyle=linestyle)
    energy_range = [-5, 6]
    efermi = [4.417301755850272, 4.920435541999894]
    shift = True
    band.multiplot(datafile, kptfile, efermi, energy_range, shift)
    fig.savefig("band.png")
