# fmt: off

'''
here the ase-abacus implementation is pasted and modified. 
Source:
https://gitlab.com/1041176461/ase-abacus/-/blob/master/ase/calculators/abacus.py

This module defines an ASE interface to ABACUS.
Created on Fri Jun  8 16:33:38 2018

ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is an open-source 
package based on density functional theory (DFT). The package utilizes both plane 
wave and numerical atomic basis sets with the usage of pseudopotentials to describe 
the interactions between nuclear ions and valence electrons. ABACUS supports LDA, 
GGA, meta-GGA, and hybrid functionals. Apart from single-point calculations, 
the package allows geometry optimizations and ab-initio molecular dynamics with 
various ensembles. The package also provides a variety of advanced functionalities 
for simulating materials, including the DFT+U, VdW corrections, and implicit solvation
model, etc. In addition, ABACUS strives to provide a general infrastructure to 
facilitate the developments and applications of novel machine-learning-assisted 
DFT methods (DeePKS, DP-GEN, DeepH, DeePTB etc.) in molecular and material simulations.

Modified on Wed Jun 20 15:00:00 2018
@author: Shen Zhen-Xiong

Modified on Wed Jun 03 23:00:00 2022
@author: Ji Yu-yang

Refactored from Sun Dec 07 21:41 2025
@author: Huang Yi-ke
'''

import os
import re
import shutil
import unittest
from pathlib import Path
from typing import Dict, Optional, List, Tuple, Set

import numpy as np
from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
    read_stdout
)
from ase.atoms import Atoms
from ase.dft.kpoints import BandPath
from ase.io import read

from abacuslite.io.generalio import (
    file_safe_backup,
    read_input,
    read_stru,
    read_kpt,
    write_input,
    write_stru,
    write_kpt
)

__LEGACYIO__ = True

class AbacusProfile(BaseProfile):
    '''AbacusProfile for interacting the ASE with ABACUS that installed in
    the practical system'''
    configvars = {'pseudo_dir', 'orbital_dir'}

    def __init__(self, 
                 command: str, 
                 pseudo_dir: Optional[str | Path] = None, 
                 orbital_dir: Optional[str | Path] = None, 
                 omp_num_threads: Optional[int] = None,
                 **kwargs):
        '''Initialize ABACUS profile.
        
        Parameters
        ----------
        command : str
            The command to run ABACUS. NOTE: there may be the case for some
            sophisticated ABACUS user they call ABACUS with command like
            `OMP_NUM_THREADS=1 mpirun -np X abacus`. Here please do not set
            the number of omp threads in `command`, instead, use `nomp=1`.
        pseudo_dir : str or Path, optional
            The directory containing pseudopotential files.
        orbital_dir : str or Path, optional
            The directory containing orbital basis files. This is only necessary
            for an ABACUS-LCAO calculation
        omp_num_threads : int, optional
            The number of omp threads to use.
        '''
        assert isinstance(command, str)
        # further validation on the command will be in the __init__ of
        # the base class
        super().__init__(command, **kwargs)
        self.pseudo_dir  = pseudo_dir
        self.orbital_dir = orbital_dir

        if omp_num_threads is not None:
            # set the number of omp threads for the present process
            assert isinstance(omp_num_threads, int)
            os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)

    @staticmethod
    def parse_version(stdout) -> str:
        # up to the ABACUS version v3.9.0.17, the run of command
        # `abacus --version` would returns the information organized
        # in the following way:
        # ABACUS version v3.9.0.17
        return re.match(r'ABACUS version (\S+)', stdout).group(1)

    def get_calculator_command(self, inputfile) -> List[str]:
        # because ABACUS run in the folder where there are INPUT files, so the
        # additional inputfile argument is not used.
        return []

    def version(self) -> str:
        '''get the abacus version information'''
        cmd_ = [*self._split_command, '--version']
        return AbacusProfile.parse_version(read_stdout(cmd_))

class AbacusTemplate(CalculatorTemplate):
    
    implemented_properties = [
        'energy', 'forces', 'stress', 'free_energy', 'magmom', 'dipole'
    ]
    _label = 'abacus'

    def __init__(self):
        super().__init__(
            'abacus',
            self.implemented_properties
        )
        self.non_convergence_ok = False
        # the redirect stdout and stderr
        self.inputname  = 'INPUT' # hard-coded
        self.outputname = f'{self._label}.out'
        self.errorname  = f'{self._label}.err'

    '''because it may be not one-to-one mapping between the property
    desired to calculate and the keywords used in the calculation,
    in the following a series of functions for mapping the property
    calculation to the keywords settings are implemented'''
    @staticmethod
    def get_energy_keywords(self) -> Dict[str, str]:
        return {}

    @staticmethod
    def get_forces_keywords(self) -> Dict[str, str]:
        return {'cal_force': '1'}
    
    @staticmethod
    def get_stress_keywords(self) -> Dict[str, str]:
        return {'cal_stress': '1'}

    @staticmethod
    def get_free_energy_keywords(self) -> Dict[str, str]:
        return {}

    @staticmethod
    def get_magmom_keywords(self) -> Dict[str, str]:
        return {'nspin': '2'}
    
    @staticmethod
    def get_dipole_keywords(self) -> Dict[str, str]:
        return {'esolver_type': 'tddft', 'out_dipole': '1'}

    def get_property_keywords(self,
                              parameters: Dict[str, str],
                              properties: List[str]) -> Dict[str, str]:
        '''Connect the relationship between the properties calculation and
        the ABACUS keywords. May be more complicated in the future, therefore
        it is better to have a seperate mapping function instead of 
        implementing in some other functions.
        
        Parameters
        ----------
        parameters : dict
            The parameters used to perform the calculation.
        properties : list of str
            The list of properties to calculate
        '''
        # update the parameters with the keywords for the properties
        # however, one should also consider that there may be the case that
        # contradictory keywords are needed. In this kind of cases, 
        # we should raise a ValueError
        param_cache_ = {}
        def counter(param_new: Dict[str, str]) -> Dict[str, str]:
            info = 'desired properties required contradictory keywords'
            for k, v in param_new.items():
                if k in param_cache_ and param_cache_[k] != v:
                    raise ValueError(f'{info}: {k}={v} (now), {param_cache_[k]} (before)')
            # if it is alright, pass through
            return param_new

        # update the parameters with the keywords for the properties
        for p in properties:
            assert p in self.implemented_properties
            parameters.update(counter(getattr(self, f'get_{p}_keywords')(parameters)))
        
        # from the parameters, get the file path
        self.suffix = parameters.get('suffix', 'ABACUS')
        self.calculation = parameters.get('calculation', 'scf')
        # with the above two, the running log file can be positioned.
        return parameters

    def write_input(self, 
                    profile: AbacusProfile, 
                    directory: Path | str,
                    atoms, 
                    parameters: Dict[str, str],
                    properties: List[str]) -> None:
        '''Write the input files for the calculation. This function connects
        the calculation in ASE language (atoms, properties, assisted by the
        parameters) to the input files of ABACUS.

        Parameters
        ----------
        profile : AbacusProfile
            The profile used to perform the calculation.
        directory : Path
            The working directory to store the input files.
        atoms : Atoms
            The atoms object to perform the calculation on. Because 
        parameters: dict
            The parameters used to perform the calculation.
        properties: list of str
            The list of properties to calculate
        '''
        # directory
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)

        # copy the `parameters` because later we will modify it
        parameters = parameters.copy()

        # STRU
        _ = file_safe_backup(directory / parameters.get('stru_file', 'STRU'))
        _ = write_stru(atoms, 
                       outdir=directory, 
                       pp_file=parameters.get('pseudopotentials'),
                       orb_file=parameters.get('basissets'),
                       fname=parameters.get('stru_file', 'STRU'))

        # KPT, if needed
        if 'kpts' in parameters:
            _ = file_safe_backup(directory / parameters.get('kpoint_file', 'KPT'))
            _ = write_kpt(parameters['kpts'], 
                          directory / parameters.get('kpoint_file', 'KPT'))
        # should this function be responsible for checking the integrity
        # of information provided by the user? There may be the case that
        # user provides incomplete information, such that the ABACUS cannot
        # run with parameters.

        # INPUT
        # after writing the KPT and STRU, delete them from the parameters
        _ = parameters.pop('kpts', None)

        _ = parameters.pop('pseudopotentials', None)
        parameters.update({'pseudo_dir': profile.pseudo_dir})

        _ = parameters.pop('basissets', None)
        parameters.update({'orbital_dir': profile.orbital_dir})
        # update the parameters respect to the properties desired
        parameters = self.get_property_keywords(parameters, properties)
        # postprocess on the parameters: convert the key and values
        # from any to string. For the case where the value is a 
        # array, convert to the string spaced by whitespace
        for k, v in parameters.items():
            # if the v is iterable, convert to the string spaced by whitespace
            if isinstance(v, (List, Tuple, Set)):
                parameters[k] = ' '.join(str(i) for i in v)
        dst = directory / self.inputname
        _ = file_safe_backup(dst)
        # remove possible key-value pairs whose value is None
        parameters = {k: v for k, v in parameters.items() if v is not None}
        _ = write_input(parameters, dst)

    def execute(self, 
                directory: Path | str, 
                profile: AbacusProfile):
        profile.run(directory=directory, 
                    inputfile=None, 
                    outputfile=self.outputname, 
                    errorfile=self.errorname)

    def read_results(self, directory) -> Dict:
        '''the function that returns the desired properties in dict'''
        read_abacus_out = lambda fn: None
        global __LEGACYIO__
        if __LEGACYIO__:
            from abacuslite.io.legacyio import read_abacus_out
        else:
            from abacuslite.io.latestio import read_abacus_out
        outdir = directory / f'OUT.{self.suffix}'
        atoms = read_abacus_out(outdir / f'running_{self.calculation}.log')[-1]
        assert atoms is not None
        return dict(atoms.calc.properties())

    def load_profile(self, cfg, **kwargs):
        return AbacusProfile.from_config(cfg, self.name, **kwargs)

class Abacus(GenericFileIOCalculator):
    def __init__(self, 
                 profile=None, 
                 directory='.', 
                 **kwargs):
        '''Construct the ABACUS calculator.

        The keyword arguments (kwargs) can be one of the ASE standard
        keywords: 'xc', 'kpts' or any of ABACUS'
        native keywords.

        Parameters
        ----------
        profile: AbacusProfile
            the interface that interacts with the ABACUS executable.
        directory: str or Path
            the working directory to store the input files.
        pseudopotentials: dict
            A mapping from the element to the pseudopotential file name,
            e.g. ``{'O': 'O_ONCV_PBE-1.0.upf', 'H': 'H.upf'}``.
        baisssets: dict, optional
            A mapping from the element to the ABACUS numerical atomic 
            orbital file name. This is necessary only when it is an
            ABACUS-LCAO (Linear-Combination-of-Atomic-Orbitals) calculation
            e.g. ``{'O': 'O_gga_10au_100Ry_2s2p1d.orb', 
            'H': 'H_gga_10au_100Ry_2s1p.orb'}``.
        kpts: dict
            The k-points sampling should be given as a dict. For there
            are many modes of k-sampling supported, the content may differ
            in cases. A `mode` key should be used to specify the ksampling,
            allowed modes are: `mp-sampling`, `line` and `point`. For 
            `mp-sampling` mode, `gamma-centered`, `nk` and `kshift` should
            present. `gamma-centered` is a boolean, `nk` and `kshift` should
            be lists of three integers. ... TBD
        inp: dict
            parameters setting in INPUT of ABACUS. NOTE: if there are settings
            on the `pseudo_dir` and `orbital_dir`, these will overwrite the
            value in the profile. If you do not expect this, please only use
            the profile, because the profile stands for interfacing with the
            ASE calculator instance with the computational environment.

        **kwargs:
            Other parameters to be passed to the ABACUS calculator.
        '''
        # not recommended :(
        profile = AbacusProfile('abacus') if profile is None else profile

        # does not support ABACUS version series v3.9.0.x
        version = profile.version()
        if re.match(r'v3\.9\.0\.\d+', version):
            global __LEGACYIO__
            __LEGACYIO__ = False
            raise NotImplementedError('ABACUS version series v3.9.0.x is not supported')

        # because ABACUS run job in folders, based on the assumption that
        # there is only one job in the folder. Therefore once there are already
        # files in the folder, will try to create a new one...(seriously?)
        inp = kwargs.pop('inp', {})

        super().__init__(
            template=AbacusTemplate(),
            profile=profile,
            parameters=kwargs | inp,
            directory=directory,
        )

    @classmethod
    def restart(cls, profile=None, directory='.', **kwargs):
        '''instantiate one ABACUS calculator from an existing job directory,
        optionally overwrite some keywords'''
        directory = Path(directory)
        inp_read = read_input(directory / 'INPUT')

        pporb_read = read_stru(directory / 'STRU')['species']
        pseudopotentials = kwargs.get(
            'pseudopotentials',
            {pporb['symbol']: pporb['pp_file'] for pporb in pporb_read}
        )
        if 'pseudopotentials' in kwargs:
            del kwargs['pseudopotentials']
        
        basissets = kwargs.get(
            'basissets',
            {pporb['symbol']: pporb.get('orb_file') for pporb in pporb_read}
        )
        if 'basissets' in kwargs:
            del kwargs['basissets']
        if all([forb is None for forb in basissets.values()]):
            basissets = {}
        assert all([forb is not None for forb in basissets.values()])

        kpts = kwargs.get('kpts', read_kpt(directory / inp_read.get('kpoint_file', 'KPT')))
        if 'kpts' in kwargs:
            del kwargs['kpts']

        inp = inp_read | kwargs.get('inp', {})
        if 'inp' in kwargs:
            del kwargs['inp']

        return cls(profile=profile, 
                   directory=directory,
                   pseudopotentials=pseudopotentials,
                   basissets=basissets,
                   kpts=kpts,
                   inp=inp,
                   **kwargs)

    def fixed_density(self,
                      kpts: BandPath | Dict[str, str | int | List[float]],
                      symmetry: str = 'off', 
                      profile=None, 
                      **kwargs) -> 'Abacus':
        '''spawn a new ABACUS calculator with fixed density, based on the present
        instance. This funcionality is mostly only useful when perform the 
        non-self-consistent calculations like band structure.
        This interface is referred from the ASE document at:
        https://ase-lib.org/gettingstarted/tut04_bulk/bulk.html#band-structure
        , however, we also note that it is from the implementation of the 
        GPAW python, not the ASE official.
        To make less development burden as possible, we use the same interface
        as the GPAW python.

        Parameters
        ----------
        kpts : BandPath | Dict[str, str | int | List[float]]
            The k-point path to be calculated. Can be either a BandPath object
            or a dictionary that contains the k-point information. For the latter
            case, see tbgen/calculators/abacus/generalio.py::write_kpt for more
            details.
        symmetry : str, optional
            The symmetry mode to be used. Default is 'off'. Now only the `off`
            mode is supported.
        profile : AbacusProfile, optional
            The profile to be used. Default is None. If None, the profile of
            the present instance will be used.
        **kwargs : dict
            Other parameters to be passed to the ABACUS calculator.
        
        Returns
        -------
        Abacus
            The new ABACUS calculator instance that can perform the nscf calculation
            tasks
        '''
        # we should overwrite the 'calculation' to 'nscf', and 'init_chg' to 'file'
        assert symmetry == 'off'
        
        kwargs.setdefault('inp', {}).update({'calculation': 'nscf',
                                             'init_chg': 'file',
                                             'symmetry': 0,
                                             'out_band': 1,
                                             'kspacing': 0.0,       # overwrite
                                             'gamma_only': False})  # overwrite

        profile = self.profile if profile is None else profile

        # get the kpoint coordinates
        if isinstance(kpts, BandPath):
            kwargs['kpts'] = {
                'mode': 'point',
                'nk': len(kpts.kpts),
                'nkinterpl': np.ones(len(kpts.kpts), dtype=int).tolist(),
                'coordinate': 'direct',
                'kpoints': kpts.kpts.tolist(),
            }
        else:
            assert isinstance(kpts, dict)
            kwargs['kpts'] = kpts
        
        # return
        return Abacus.restart(profile=profile, 
                              directory=self.directory,
                              **kwargs)

    def band_structure(self, efermi=None):
        '''get the band structure from ABACUS. 
        (now not only GPAW can calculate the band structure ;) )'''
        from ase.spectrum.band_structure import get_band_structure
        return get_band_structure(calc=self, reference=efermi)

class TestAbacusCalculator(unittest.TestCase):

    here = Path(__file__).parent
    testfiles = here.parent.parent / 'testfiles'

    @unittest.skip('This is an example to perform single-point calculation')
    def test_get_potential_energy(self):
        from ase.build.bulk import bulk
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 8 abacus',
            pseudo_dir=self.testfiles / 'pporb',
            orbital_dir=self.testfiles / 'pporb',
            omp_num_threads=1
        )
        calculator = Abacus(aprof,
                            directory=self.here / 'test_get_potential_energy',
                            pseudopotentials={'Si': 'Si.upf'},
                            basissets={'Si': 'Si_gga_7au_100Ry_2s2p1d.orb'},
                            inp={'calculation': 'scf',
                                 'basis_type': 'lcao',
                                 'ks_solver': 'genelpa',
                                 'ecutwfc': 40,
                                 'symmetry': 1,
                                 'nspin': 1,
                                 'cal_force': 1,
                                 'cal_stress': 1},
                            kpts={'mode': 'mp-sampling',
                                  'gamma-centered': True,
                                  'nk': (3, 3, 3),
                                  'kshift': (0, 0, 0)})
        silicon.calc = calculator
        e = silicon.get_potential_energy()
        self.assertAlmostEqual(e, -229.5298358034765727, places=8)
        # print(calculator.results)
        # {
        #   'nspins': 1, 
        #   'nkpts': 4, 
        #   'nbands': 14, 
        #   'eigenvalues': array([[[-5.64733 ,  6.34567 ,  6.34567 ,  6.34567 ,  8.87785 ,
        #                            8.87785 ,  8.87785 ,  9.71303 , 15.083   , 15.083   ,
        #                            22.9239  , 22.9239  , 22.9239  , 31.3422  ],
        #                          [-4.27085 ,  0.981027,  5.32824 ,  5.32824 ,  8.14638 ,
        #                            9.93948 ,  9.93948 , 14.7114  , 16.2124  , 16.2124  ,
        #                           21.9002  , 22.2073  , 22.213   , 22.213   ],
        #                          [-3.74469 ,  1.31154 ,  3.90001 ,  3.90001 ,  7.05826 ,
        #                            8.82102 , 13.6787  , 13.6787  , 17.6838  , 18.2193  ,
        #                           20.5768  , 20.5768  , 21.9485  , 25.0883  ],
        #                          [-2.38241 , -0.457849,  1.8921  ,  4.20947 ,  7.89527 ,
        #                           12.7831  , 12.8361  , 13.857   , 14.1422  , 15.7449  ,
        #                           22.4027  , 23.5822  , 24.1805  , 24.2308  ]]]), 
        #   'occupations': array([[[0.0740741 , 0.0733693 , 0.0733693 , 0.0733693 , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        ],
        #                          [0.592593  , 0.592593  , 0.592593  , 0.592593  , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        ],
        #                          [0.444444  , 0.444444  , 0.444444  , 0.444444  , 0.00211439,
        #                           0.        , 0.        , 0.        , 0.        , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        ],
        #                          [0.888889  , 0.888889  , 0.888889  , 0.888889  , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        , 0.        ,
        #                           0.        , 0.        , 0.        , 0.        ]]]), 
        #   'fermi_level': 6.6840689668, 
        #   'kpoint_weights': array([0.037 , 0.2963, 0.2222, 0.4444]), 
        #   'ibz_kpoints': array([[ 0.        ,  0.        ,  0.        ],
        #                         [ 0.33333333,  0.33333333,  0.33333333],
        #                         [ 0.33333333,  0.33333333,  0.        ],
        #                         [-0.33333333,  0.33333333,  0.33333333]]), 
        #   'energy': -229.5298358035, 
        #   'free_energy': -229.5298358035, 
        #   'natoms': 2, 
        #   'forces': array([[0., 0., 0.],
        #                    [0., 0., 0.]]), 
        #   'stress': array([44.08253252, 44.08253252, 44.08253252,
        #                     0.        ,  0.        , -0.        ]), 
        #   'magmoms': array([0., 0.])
        # }
        
        # remove the jobdir
        shutil.rmtree(self.here / 'test_get_potential_energy')

    @unittest.skip('This is an example to calculate the DOS')
    def test_dos(self):
        '''try if possible to integrate the ABACUS-ASE Plugin to calculate the density of states
        reference: https://ase-lib.org/gettingstarted/tut04_bulk/bulk.html#density-of-states
        '''
        from ase.dft.dos import DOS
        from ase.build.bulk import bulk
        import matplotlib.pyplot as plt

        # first perform the SCF calculation
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 8 abacus',
            pseudo_dir=self.testfiles / 'pporb',
            orbital_dir=self.testfiles / 'pporb',
            omp_num_threads=1
        )
        calculator = Abacus(aprof,
                            directory=self.here / 'test_get_potential_energy',
                            pseudopotentials={'Si': 'Si.upf'},
                            basissets={'Si': 'Si_gga_7au_100Ry_2s2p1d.orb'},
                            inp={'calculation': 'scf',
                                 'basis_type': 'lcao',
                                 'ks_solver': 'genelpa',
                                 'ecutwfc': 40,
                                 'symmetry': 1,
                                 'nspin': 1,
                                 'cal_force': 1,
                                 'cal_stress': 1},
                            kpts={'mode': 'mp-sampling',
                                  'gamma-centered': True,
                                  'nk': (3, 3, 3),
                                  'kshift': (0, 0, 0)})
        silicon.calc = calculator
        e = silicon.get_potential_energy()

        # then analyze the DOS
        dos = DOS(calculator, width=0.1)
        energies = dos.get_energies()
        weights = dos.get_dos()
        plt.plot(energies, weights)
        plt.show()

    @unittest.skip('This is an example to calculate the bandstructure')
    def test_bandstructure(self):
        '''try if possible to integrate the ABACUS-ASE Plugin to calculate the band structure
        reference: https://ase-lib.org/gettingstarted/tut04_bulk/bulk.html#band-structure
        '''
        import seekpath
        from abacuslite.utils.ksampling import merge_ksgm, make_kstring

        # first perform the SCF calculation
        silicon: Atoms = read(self.testfiles / 'Si_mp-149_computed.cif')
        aprof = AbacusProfile(
            command='mpirun -np 8 abacus',
            pseudo_dir=self.testfiles / 'pporb',
            orbital_dir=self.testfiles / 'pporb',
            omp_num_threads=1
        )
        calculator = Abacus(aprof,
                            directory=self.here / 'test_bandstructure',
                            pseudopotentials={'Si': 'Si.upf'},
                            basissets={'Si': 'Si_gga_7au_100Ry_2s2p1d.orb'},
                            inp={'calculation': 'scf',
                                 'basis_type': 'lcao',
                                 'ks_solver': 'genelpa',
                                 'ecutwfc': 40,
                                 'symmetry': 1,
                                 'nspin': 1,
                                 'cal_force': 1,
                                 'cal_stress': 1,
                                 'out_chg': '1 8',
                                 'kspacing': 0.1}, # for the band structure calculation
                            )
        silicon.calc = calculator
        _ = silicon.get_potential_energy()

        # then the band structure non-self-consistent calculation
        kpathseen = seekpath.get_path(
            structure=(np.array(silicon.get_cell()), 
                       silicon.get_scaled_positions(), 
                       silicon.get_atomic_numbers()),
            with_time_reversal=True
        )

        # seems the following lines should be moved elsewhere...? but where?
        kpathstr, is_brkpt = merge_ksgm(kpathseen['path'])
        fklblfilter = lambda lbl: lbl if lbl != 'GAMMA' else 'G'
        kpathstr = make_kstring([fklblfilter(lbl) for lbl in kpathstr], is_brkpt)
        kpathstr = ''.join([k if k != 'GAMMA' else 'G' for k in kpathstr])
        # seekpath use 'GAMMA' as the gamma point symbol, while the ASE use 'G'
        kspecial = {k: v for k, v in kpathseen['point_coords'].items() if k != 'GAMMA'}
        kspecial['G'] = [0.0, 0.0, 0.0]

        bandpath = silicon.cell.bandpath(path=kpathstr,
                                         npoints=50,
                                         special_points=kspecial)
        scfcalc: Abacus = silicon.calc
        bscalc = scfcalc.fixed_density(bandpath)
        silicon.calc = bscalc
        _ = silicon.get_potential_energy() # let's see what happened here

        # remove the files...
        shutil.rmtree(self.here / 'test_bandstructure')

        eigenvalues = bscalc.results['eigenvalues']
        self.assertIsInstance(eigenvalues, np.ndarray)
        self.assertEqual(eigenvalues.shape, (1, 50, 14)) # 50 kpoints, 14 bands

        # postprocessing the band structure
        # bs = bscalc.band_structure()
        # bs.write('bs.json')
        # then in the terminal, type the following command
        # ```
        # ase band-structure bs.json -r -20 20
        # ```
        # to plot the band structure

    @unittest.skip('Not completed yet')
    def test_restart(self):
        aprof = AbacusProfile(
            command='mpirun -np 8 abacus',
            pseudo_dir=self.testfiles / 'pporb',
            orbital_dir=self.testfiles / 'pporb',
            omp_num_threads=1
        )
        calc = Abacus.restart(aprof,
                              directory=self.testfiles / 'sih-1d')
        bscalc = calc.fixed_density()

if __name__ == '__main__':
    unittest.main()