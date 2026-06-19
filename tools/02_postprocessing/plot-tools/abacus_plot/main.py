'''
Date: 2022-01-03 15:03:54
LastEditors: jiyuyang
LastEditTime: 2022-01-03 17:10:52
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import argparse
import matplotlib.pyplot as plt

from abacus_plot.band import Band, PBand
from abacus_plot.dos import TDOS, PDOS
from abacus_plot.utils import read_json, key2int


class Show:
    """Show auto-test information"""

    @classmethod
    def show_cmdline(cls, args):
        if args.band and not args.projected and not args.out:
            text = read_json(args.file)
            bandfile = text["bandfile"]
            kptfile = text["kptfile"]
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            shift = text.pop("shift", False)
            figsize = text.pop("figsize", (12, 10))
            fig, ax = plt.subplots(figsize=figsize)
            band = Band(bandfile, kptfile)
            band.plot(fig=fig, ax=ax, efermi=efermi, energy_range=energy_range, shift=shift, **text)
            dpi = text.pop("dpi", 400)
            fig.savefig(text.pop("outfile", "band.png"), dpi=dpi)

        if args.band and args.projected:
            text = read_json(args.file)
            bandfile = text["bandfile"]
            kptfile = text["kptfile"]
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            shift = text.pop("shift", False)
            figsize = text.pop("figsize", (12, 10))
            fig, ax = plt.subplots(figsize=figsize)
            index = key2int(text.pop("index", None))    # keys of json must be dict, I convert them to int for `index` and `atom_index`
            atom_index = key2int(text.pop("atom_index", None))
            species = text.pop("species", None)
            outdir = text.pop("outdir", './')
            cmapname = text.pop("cmapname", 'jet')
            pband = PBand(bandfile, kptfile)
            pband.plot(fig=fig, ax=ax, index=index, atom_index=atom_index, species=species, efermi=efermi, energy_range=energy_range, shift=shift, outdir=outdir, cmapname=cmapname, **text)

        if args.band and args.out:
            text = read_json(args.file)
            bandfile = text["bandfile"]
            kptfile = text["kptfile"]
            index = key2int(text.pop("index", None))
            atom_index = key2int(text.pop("atom_index", None))
            species = text.pop("species", None)
            outdir = text.pop("outdir", './')
            pdos = PBand(bandfile, kptfile)
            pdos.write(index=index, atom_index=atom_index, species=species, outdir=outdir)

        if args.dos and not args.projected and not args.out:
            text = read_json(args.file)
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

        if args.dos and args.projected:
            text = read_json(args.file)
            pdosfile = text.pop("pdosfile", '')
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", [])
            dos_range = text.pop("dos_range", [])
            shift = text.pop("shift", False)
            index = key2int(text.pop("index", None))
            atom_index = key2int(text.pop("atom_index", None))
            species = text.pop("species", None)
            figsize = text.pop("figsize", (12, 10))
            if index:
                fig, ax = plt.subplots(
                    len(index), 1, sharex=True, figsize=figsize)
            if atom_index:
                fig, ax = plt.subplots(
                    len(atom_index), 1, sharex=True, figsize=figsize)
            if species:
                fig, ax = plt.subplots(
                    len(species), 1, sharex=True, figsize=figsize)
            prec = text.pop("prec", 0.01)
            pdos = PDOS(pdosfile)
            pdos.plot(fig=fig, ax=ax, index=index, atom_index=atom_index, species=species, efermi=efermi, shift=shift,
                      energy_range=energy_range, dos_range=dos_range, prec=prec, **text)
            dpi = text.pop("dpi", 400)
            fig.savefig(text.pop("pdosfig", "pdos.png"), dpi=dpi)

        if args.dos and args.out:
            text = read_json(args.file)
            pdosfile = text.pop("pdosfile", '')
            index = key2int(text.pop("index", None))
            atom_index = key2int(text.pop("atom_index", None))
            species = text.pop("species", None)
            outdir = text.pop("outdir", './')
            pdos = PDOS(pdosfile)
            pdos.write(index=index, atom_index=atom_index, species=species, outdir=outdir)

def main():
    parser = argparse.ArgumentParser(
        prog='abacus-plot', description='Plotting tools for ABACUS')

    # Show
    parser.add_argument('-f', '--file', dest='file', type=str, nargs='?', const='config.json',
                        default='config.json', help='profile with format json')
    parser.add_argument('-b', '--band', dest='band', nargs='?', const=1, type=int,
                        default=None, help='plot band structure and show band information.')
    parser.add_argument('-d', '--dos', dest='dos', nargs='?', const=1, type=int,
                        default=None, help='plot density of state(DOS).')
    parser.add_argument('-p', '--projected', dest='projected', nargs='?', const=1, type=int,
                        default=None, help='plot projected band structure or partial density of state(PDOS), should be used with `-b` or `-d`.')
    parser.add_argument('-o', '--out_parsed_data', dest='out', nargs='?', const=1, type=int,
                        default=None, help='output projected band structure or partial density of state(PDOS) to files, should be used with `-b` or `-d`.')
    parser.set_defaults(func=Show().show_cmdline)

    args = parser.parse_args()
    args.func(args)
