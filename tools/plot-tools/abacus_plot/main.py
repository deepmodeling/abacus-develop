'''
Date: 2022-01-03 15:03:54
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:10:52
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import argparse
from os import PathLike
import matplotlib.pyplot as plt

from abacus_plot.band import BandPlot
from abacus_plot.dos import TDOS, PDOS
from abacus_plot.utils import read_json


class Show:
    """Show auto-test information"""

    @classmethod
    def show_cmdline(cls, args):
        if args.band:
            text = read_json(args.band)
            datafile = text["bandfile"]
            kptfile = text["kptfile"]
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            shift = text.pop("shift", False)
            figsize = text.pop("figsize", (12, 10))
            fig, ax = plt.subplots(figsize=figsize)
            band = BandPlot(fig, ax)
            if isinstance(datafile, (str, PathLike)):
                band.singleplot(datafile, kptfile, efermi, energy_range, shift)
            elif isinstance(datafile, (list, tuple)):
                band.multiplot(datafile, kptfile, efermi, energy_range, shift)
            dpi = text.pop("dpi", 400)
            fig.savefig(text.pop("outfile", "band.png"), dpi=dpi)

        if args.tdos:
            text = read_json(args.tdos)
            tdosfile = text.pop("tdosfile", '')
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            dos_range = text.pop("dos_range", [])
            shift = text.pop("shift", False)
            figsize = text.pop("figsize", (12, 10))
            fig, ax = plt.subplots(figsize=figsize)
            prec = text.pop("prec", 0.01)
            tdos = TDOS(tdosfile)
            tdos.plot(fig, ax, efermi=efermi, shift=shift,
                      energy_range=energy_range, dos_range=dos_range, prec=prec, **text)
            dpi = text.pop("dpi", 400)
            fig.savefig(text.pop("tdosfig", "tdos.png"), dpi=dpi)

        if args.pdos:
            text = read_json(args.pdos)
            pdosfile = text.pop("pdosfile", '')
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            dos_range = text.pop("dos_range", [])
            shift = text.pop("shift", False)
            species = text.pop("species", [])
            figsize = text.pop("figsize", (12, 10))
            fig, ax = plt.subplots(
                len(species), 1, sharex=True, figsize=figsize)
            prec = text.pop("prec", 0.01)
            pdos = PDOS(pdosfile)
            pdos.plot(fig, ax, species, efermi=efermi, shift=shift,
                      energy_range=energy_range, dos_range=dos_range, prec=prec, **text)
            dpi = text.pop("dpi", 400)
            fig.savefig(text.pop("pdosfig", "pdos.png"), dpi=dpi)

        if args.out_pdos:
            text = read_json(args.out_pdos)
            pdosfile = text.pop("pdosfile", '')
            species = text.pop("species", [])
            outdir = text.pop("outdir", './')
            pdos = PDOS(pdosfile)
            pdos.write(species, outdir)

def main():
    parser = argparse.ArgumentParser(
        prog='abacus-plot', description='Plotting tools for ABACUS')

    # Show
    parser.add_argument('-b', '--band', dest='band', type=str,
                        default=None, help='plot band structure and show band information.')
    parser.add_argument('-t', '--tdos', dest='tdos', type=str,
                        default=None, help='plot total density of state(TDOS).')
    parser.add_argument('-p', '--pdos', dest='pdos', type=str,
                        default=None, help='plot partial density of state(PDOS).')
    parser.add_argument('-o', '--out_pdos', dest='out_pdos', type=str,
                        default=None, help='output partial density of state(PDOS).')
    parser.set_defaults(func=Show().show_cmdline)

    args = parser.parse_args()
    args.func(args)
