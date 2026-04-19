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
import tempfile
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

        # fix: inconsistent atoms order may induce bugs, here a list
        # is kept to swap the order of atoms
        self.atomorder  = None

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
                    atoms: Atoms, 
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
        # reorder the atoms according to the alphabet. Keep the reverse map
        # so that we will recover the order in function read_results()
        ind = sorted(range(len(atoms)), key=lambda i: atoms[i].symbol)
        self.atomorder = sorted(range(len(atoms)), key=lambda i: ind[i]) # revmap
        # then we write
        _ = write_stru(atoms[ind], 
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

        # FIXME: only support the ksdft esolver_type presently
        if parameters.get('esolver_type', 'ksdft') != 'ksdft':
            raise NotImplementedError(
                'ABACUS Lite only supports the ksdft esolver_type presently, '
                'which means the ABACUS should always be used as a DFT '
                'calculator. For other forcefields that ABACUS supports '
                'such as the LJ, DP, etc., please either use the ABACUS '
                'directly, or the implementation of interfaces to ASE '
                'directly.'
            )

        # write the INPUT file to the target directory
        _ = write_input(parameters, dst)

    def execute(self, 
                directory: Path | str, 
                profile: AbacusProfile):
        '''Execute the ABACUS Lite calculation.

        Parameters
        ----------
        directory : Path or str
            The working directory to store the input files.
        profile : AbacusProfile
            The profile used to perform the calculation.

        Raises
        ------
        SubprocessError
            If the ABACUS Lite calculation fails.
        '''
        from subprocess import SubprocessError
        try:
            profile.run(directory=directory, 
                        inputfile=None, 
                        outputfile=self.outputname, 
                        errorfile=self.errorname)
        except SubprocessError:
            message = ['ABACUS Lite calculation failed']
            with open(directory / self.outputname, 'r') as f:
                message.append(f.read())
            with open(directory / self.errorname, 'r') as f:
                message.append(f.read())
            raise SubprocessError('\n'.join(message))

    def read_results(self, directory) -> Dict:
        '''the function that returns the desired properties in dict'''
        read_abacus_out = lambda fn: None
        global __LEGACYIO__
        if __LEGACYIO__:
            from abacuslite.io.legacyio import read_abacus_out
        else:
            from abacuslite.io.latestio import read_abacus_out

        outdir = directory / f'OUT.{self.suffix}'
        # only the last frame
        atoms: Optional[Atoms] = read_abacus_out(
            outdir / f'running_{self.calculation}.log',
            sort_atoms_with=self.atomorder)[-1]
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

        # does not support ABACUS version series v3.9.0.x and v3.11.0-beta.x
        version = profile.version()
        if re.match(r'v3\.9\.0\.\d+', version) or re.match(r'v3\.11\.0-beta\.\d+', version):
            global __LEGACYIO__
            __LEGACYIO__ = False

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
    pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

    def test_calculator_results(self):
        from ase.build.bulk import bulk
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=self.pporb,
            orbital_dir=self.pporb,
            omp_num_threads=1
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            calculator = Abacus(aprof,
                                directory=tmpdir,
                                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                                basissets={'Si': 'Si_gga_6au_100Ry_2s2p1d.orb'},
                                inp={'calculation': 'scf',
                                    'basis_type': 'lcao',
                                    'ks_solver': 'genelpa',
                                    'ecutwfc': 40,
                                    'symmetry': 1,
                                    'nspin': 1,
                                    'gamma_only': True,
                                    'cal_force': 1,
                                    'cal_stress': 1})
            silicon.calc = calculator
            e = silicon.get_potential_energy()
        
        # check!
        self.assertAlmostEqual(e, -194.953053309)
        self.assertIsNotNone(calculator.results)
        self.assertIsInstance(calculator.results, dict)
        for k in ['nspins', 'nkpts', 'nbands', 'eigenvalues', 'occupations',
                  'fermi_level', 'kpoint_weights', 'ibz_kpoints', 'energy', 
                  'free_energy', 'natoms', 'forces', 'stress', 'magmoms']:
            self.assertIn(k, calculator.results)
        self.assertEqual(calculator.results['nspins'], 1)
        self.assertEqual(calculator.results['nkpts'], 1)
        self.assertEqual(calculator.results['nbands'], 14)
        self.assertEqual(calculator.results['energy'], e)
        self.assertEqual(calculator.results['free_energy'], e)
        self.assertEqual(calculator.results['natoms'], 2)
        
        for k in ['eigenvalues', 'occupations', 'ibz_kpoints', 'forces', 'stress', 'magmoms']:
            self.assertIsInstance(calculator.results[k], np.ndarray)

        self.assertEqual(calculator.results['eigenvalues'].shape, (1, 1, 14))
        ekb = [-4.82194,  7.62727,  7.62727,  7.62737, 10.2436 , 10.2436 ,
                10.2436 , 10.9884 , 16.057  , 16.057  , 23.8353 , 25.421  ,
                25.421  , 25.4212 ]
        self.assertTrue(np.allclose(calculator.results['eigenvalues'][0, 0, :], np.array(ekb)))

        self.assertEqual(calculator.results['occupations'].shape, (1, 1, 14))
        occ = [2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.assertTrue(np.allclose(calculator.results['occupations'][0, 0, :], np.array(occ)))
        
        self.assertEqual(calculator.results['ibz_kpoints'].shape, (1, 3))
        self.assertTrue(np.allclose(calculator.results['ibz_kpoints'][0, :], np.array([0,0,0])))

        self.assertEqual(calculator.results['forces'].shape, (2, 3))
        self.assertTrue(np.allclose(calculator.results['forces'], np.zeros((2, 3))))

        self.assertEqual(calculator.results['stress'].shape, (6,))
        stress = [-0.19327923, -0.19327923, -0.19327923, -0.        ,  0.        ,   0.        ]
        self.assertTrue(np.allclose(calculator.results['stress'], np.array(stress)))
        
        self.assertEqual(calculator.results['magmoms'].shape, (2,))
        self.assertTrue(np.allclose(calculator.results['magmoms'], np.zeros(2)))

    def test_restart(self):
        from ase.build.bulk import bulk
        silicon = bulk('Si', crystalstructure='diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=self.pporb,
            orbital_dir=self.pporb,
            omp_num_threads=1
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            calculator = Abacus(aprof,
                                directory=tmpdir,
                                pseudopotentials={'Si': 'Si_ONCV_PBE-1.0.upf'},
                                basissets={'Si': 'Si_gga_6au_100Ry_2s2p1d.orb'},
                                inp={'calculation': 'scf',
                                    'basis_type': 'lcao',
                                    'ks_solver': 'genelpa',
                                    'ecutwfc': 40,
                                    'symmetry': 1,
                                    'nspin': 1,
                                    'gamma_only': True,
                                    'cal_force': 1,
                                    'cal_stress': 1})
            silicon.calc = calculator
            e = silicon.get_potential_energy()
        
            # restart
            silicon.calc = Abacus.restart(aprof, directory=tmpdir)
            e2 = silicon.get_potential_energy()
            self.assertAlmostEqual(e2, e)

if __name__ == '__main__':
    unittest.main()