#!/usr/bin/env python3
"""Positive and corruption controls for the full PW eigenvector comparison."""
import cmath
import copy
import json
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

_spec = importlib.util.spec_from_file_location(
    "compare_pw_wfc", Path(__file__).with_name("compare_pw_wfc.py"))
comparator = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(comparator)
compare_directories = comparator.compare_directories
compare_wavefunctions = comparator.compare_wavefunctions
read_wavefunction = comparator.read_wavefunction


def fixture(index=1, count=2, npol=1):
    vectors = [(0,0,0),(0,0,1),(0,1,0),(1,0,0)]
    size = len(vectors) * npol
    bands = [[cmath.exp(2j*math.pi*band*i/size)/math.sqrt(size) for i in range(size)]
             for band in range(3)]
    return {'header': {'Kpoint': index, 'nKpoint': count, 'kv.x': 0.1*(index-1),
            'kv.y': 0., 'kv.z': 0., 'weight': 0.5, 'ngtot': 4, 'nband': 3,
            'ecut': 20., 'lat0': 10., '2pi/lat0': 0.628318},
            'lattice': [(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)], 'vectors': vectors,
            'npol': npol, 'bands': bands}


def text(wave):
    header = wave['header']
    keys = 'Kpoint nKpoint kv.x kv.y kv.z weight ngtot nband ecut lat0 2pi/lat0'.split()
    lines = [' '.join(keys), ' '.join(str(header[key]) for key in keys),
             '<Reciprocal Lattice Vector>']
    lines += [' '.join(str(x) for x in row) for row in wave['lattice']]
    lines += ['<Reciprocal Lattice Vector>', '<G vectors>']
    lines += [' '.join(str(x) for x in row) for row in wave['vectors']]
    lines += ['<G vectors>']
    for index, band in enumerate(wave['bands'],1):
        marker = f'< Band {index} >'
        lines.append(marker)
        lines.extend(f'{z.real:.16e} {z.imag:.16e}' for z in band)
        lines.append(marker)
    return '\n'.join(lines)+'\n'


class WavefunctionComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.ref = self.root/'reference'
        self.cal = self.root/'actual'
        self.ref.mkdir()
        self.cal.mkdir()
        for index in [1,2]:
            name = f'wfk{index}_pw.txt'
            for directory in [self.ref,self.cal]:
                (directory/name).write_text(text(fixture(index)))
        bands = '1 0 -2 -0.5 2\n2 0.1 -1.9 -0.4 2.1\n'
        for directory in [self.ref,self.cal]:
            (directory/'band.txt').write_text(bands)

    def tearDown(self):
        self.temporary.cleanup()

    def compare(self):
        return compare_directories(self.ref,self.cal,self.ref/'band.txt',
                                   self.cal/'band.txt',1e-7,1e-8)

    def write_actual(self,wave):
        name=f'wfk{wave["header"]["Kpoint"]}_pw.txt'
        (self.cal/name).write_text(text(wave))

    def test_identical_full_files(self):
        result=self.compare()
        self.assertEqual(len(result['kpoints']),2)
        self.assertEqual(result['max_energy_error_eV'],0)

    def test_independent_global_complex_phase_each_band_and_kpoint(self):
        for index in [1,2]:
            wave=fixture(index)
            for ib,band in enumerate(wave['bands']):
                phase=cmath.exp(1j*(ib+1)*index*0.37)
                wave['bands'][ib]=[phase*z for z in band]
            self.write_actual(wave)
        self.compare()

    def test_valid_phase_changes_the_legacy_maximum_statistic(self):
        left=fixture(1);right=copy.deepcopy(left)
        right['bands'][0]=[z*cmath.exp(0.37j) for z in right['bands'][0]]
        def legacy_maximum(wave):
            return round(max(max(abs(z.real),abs(z.imag)) for z in wave['bands'][0]),4)
        self.assertNotEqual(legacy_maximum(left),legacy_maximum(right))
        compare_wavefunctions(left,right,1e-7)

    def test_g_vector_order_matches_coefficients(self):
        wave=fixture(2)
        order=[2,0,3,1]
        wave['vectors']=[wave['vectors'][i] for i in order]
        wave['bands']=[[band[i] for i in order] for band in wave['bands']]
        self.write_actual(wave)
        self.compare()

    def test_relative_phase_corruption_preserving_old_maximum_is_rejected(self):
        wave=fixture(1)
        old_max=max(max(abs(z.real),abs(z.imag)) for z in wave['bands'][0])
        wave['bands'][0][0]*=cmath.exp(0.1j)
        self.assertEqual(old_max,max(max(abs(z.real),abs(z.imag)) for z in wave['bands'][0]))
        self.write_actual(wave)
        with self.assertRaises(ValueError):self.compare()

    def test_nondegenerate_band_rotation_is_rejected(self):
        wave=fixture(1)
        a,b=copy.deepcopy(wave['bands'][:2])
        angle=0.02
        wave['bands'][0]=[math.cos(angle)*x+math.sin(angle)*y for x,y in zip(a,b)]
        wave['bands'][1]=[-math.sin(angle)*x+math.cos(angle)*y for x,y in zip(a,b)]
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'coefficient error'):self.compare()

    def test_corruption_at_second_kpoint_is_rejected(self):
        wave=fixture(2);wave['bands'][2][3]+=0.001
        self.write_actual(wave)
        with self.assertRaises(ValueError):self.compare()

    def test_scale_is_not_normalized_away(self):
        wave=fixture(1);wave['bands'][0]=[1.001*z for z in wave['bands'][0]]
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'normalization'):self.compare()

    def test_matching_invalid_normalization_is_rejected(self):
        wave=fixture(1);wave['bands'][0]=[2*z for z in wave['bands'][0]]
        with self.assertRaisesRegex(ValueError,'normalization'):
            compare_wavefunctions(wave,wave,1e-7)

    def test_band_energy_tolerance_is_enforced(self):
        (self.cal/'band.txt').write_text('1 0 -2 -0.5 2\n2 0.1 -1.9 -0.4 2.10000002\n')
        with self.assertRaisesRegex(ValueError,'band energy error'):self.compare()

    def test_degenerate_spectrum_does_not_silently_pass(self):
        for directory in [self.ref,self.cal]:
            (directory/'band.txt').write_text('1 0 -2 -2 2\n2 0.1 -1.9 -0.4 2.1\n')
        with self.assertRaisesRegex(ValueError,'separate subspace comparison'):self.compare()

    def test_truncated_coefficients(self):
        wave=fixture(1);wave['bands'][0].pop()
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'coefficient count'):self.compare()

    def test_duplicate_g_vector(self):
        wave=fixture(1);wave['vectors'][1]=wave['vectors'][0]
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'duplicate G'):self.compare()

    def test_g_vector_corruption(self):
        wave=fixture(1);wave['vectors'][0]=(4,4,4)
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'vectors mismatch'):self.compare()

    def test_nonfinite_coefficient(self):
        wave=fixture(1);wave['bands'][0][0]=complex(float('nan'),0)
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'nonfinite'):self.compare()

    def test_nonfinite_metadata(self):
        wave=fixture(1);wave['header']['weight']=float('inf')
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'nonfinite'):self.compare()

    def test_metadata_mismatch(self):
        wave=fixture(1);wave['header']['ecut']=21
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'header mismatch'):self.compare()

    def test_missing_band(self):
        wave=fixture(1);wave['bands'].pop()
        self.write_actual(wave)
        with self.assertRaisesRegex(ValueError,'truncated'):self.compare()

    def test_unexpected_trailing_content(self):
        path=self.cal/'wfk1_pw.txt';path.write_text(path.read_text()+'unexpected 1 2\n')
        with self.assertRaisesRegex(ValueError,'trailing'):self.compare()

    def test_missing_kpoint(self):
        (self.cal/'wfk2_pw.txt').unlink()
        with self.assertRaisesRegex(ValueError,'incomplete k-point'):self.compare()

    def test_duplicate_kpoint(self):
        (self.cal/'wfk9_pw.txt').write_text(text(fixture(1)))
        with self.assertRaisesRegex(ValueError,'duplicate k-point'):self.compare()

    def test_empty_directory(self):
        empty=self.root/'empty';empty.mkdir()
        with self.assertRaisesRegex(ValueError,'no PW wavefunctions'):
            compare_directories(empty,self.cal,empty/'band.txt',self.cal/'band.txt',1e-7,1e-8)

    def test_all_floating_tolerances_must_be_finite_positive(self):
        for value in [0,-1,float('nan'),float('inf')]:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    compare_directories(self.ref,self.cal,self.ref/'band.txt',self.cal/'band.txt',value,1e-8)
                with self.assertRaises(ValueError):
                    compare_directories(self.ref,self.cal,self.ref/'band.txt',self.cal/'band.txt',1e-7,value)

    def test_spinor_phase_is_shared_by_both_components(self):
        left=fixture(1,npol=2);right=copy.deepcopy(left)
        right['bands']=[[z*cmath.exp(0.4j) for z in band] for band in right['bands']]
        compare_wavefunctions(left,right,1e-7)
        right['bands'][0][4:]=[z*cmath.exp(0.3j) for z in right['bands'][0][4:]]
        with self.assertRaises(ValueError):compare_wavefunctions(left,right,1e-7)

    def test_precision_roundtrip_fixture_contains_nontrivial_doubles(self):
        wave=fixture(1,npol=2)
        for band in wave['bands']:
            for i,z in enumerate(band):band[i]=z*cmath.exp(0.6789012345678912j)
        p=self.root/'spinor.txt';p.write_text(text(wave));read=read_wavefunction(p)
        # G sorting is independent of storage order; compare in the parsed order.
        order=sorted(range(4),key=wave['vectors'].__getitem__)
        expected=[[band[s*4+i] for s in range(2) for i in order] for band in wave['bands']]
        self.assertEqual(expected,read['bands'])

    def test_cli_success_and_nonzero_rejection(self):
        command=[sys.executable,str(Path(__file__).with_name('compare_pw_wfc.py')),str(self.ref),
                 str(self.cal),str(self.ref/'band.txt'),str(self.cal/'band.txt')]
        result=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        lines=result.stdout.splitlines()
        expected={f'Read {path.resolve()}' for directory in (self.ref,self.cal)
                  for path in directory.iterdir()}
        self.assertEqual(set(lines[:-1]),expected)
        self.assertEqual(len(json.loads(lines[-1])['kpoints']),2)
        wave=fixture(2);wave['bands'][1][2]+=0.01;self.write_actual(wave)
        result=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(result.returncode,1)
        self.assertIn('comparison failed',result.stderr)

    def run_collector(self, full_reference=True):
        directory=self.root/'collector'
        directory.mkdir()
        import shutil
        shutil.copytree(self.cal,directory/'OUT.autotest')
        if full_reference:
            shutil.copytree(self.ref,directory/'wfc_pw.ref')
            shutil.copy2(self.ref/'band.txt',directory/'band.txt.ref')
            (directory/'result.ref').write_text('CompareWfc_pass 0\n')
        tools=Path(__file__).resolve().parent
        script='source "$1/props_cube.sh"; PROPS_TOOLS_DIR="$1"; props_result_file=result.out; has_wfc_pw=1; nspin=1; run_cube_props'
        result=subprocess.run(['bash','-e','-c',script,'collector',str(tools)],cwd=directory,capture_output=True,text=True)
        return directory,result

    def test_collector_full_reference_passes(self):
        directory,result=self.run_collector()
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual((directory/'result.out').read_text(),'CompareWfc_pass 0\n')

    def test_collector_corruption_returns_failure(self):
        wave=fixture(2);wave['bands'][1][2]+=0.01;self.write_actual(wave)
        directory,result=self.run_collector()
        self.assertNotEqual(result.returncode,0)
        self.assertEqual((directory/'result.out').read_text(),'CompareWfc_pass 1\n')

    def test_collector_missing_full_reference_cannot_fall_back(self):
        directory,result=self.run_collector()
        self.assertEqual(result.returncode,0,result.stderr)
        import shutil
        shutil.rmtree(directory/'wfc_pw.ref')
        (directory/'result.out').unlink()
        tools=Path(__file__).resolve().parent
        script='source "$1/props_cube.sh"; PROPS_TOOLS_DIR="$1"; props_result_file=result.out; has_wfc_pw=1; nspin=1; run_cube_props'
        result=subprocess.run(['bash','-e','-c',script,'collector',str(tools)],cwd=directory,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual((directory/'result.out').read_text(),'CompareWfc_pass 1\n')

    def test_collector_required_reference_with_disabled_or_absent_output(self):
        directory,result=self.run_collector()
        self.assertEqual(result.returncode,0,result.stderr)
        import shutil
        shutil.rmtree(directory/'wfc_pw.ref')
        tools=Path(__file__).resolve().parent
        for flag in ['0','']:
            with self.subTest(output_flag=flag):
                (directory/'result.out').unlink()
                script='source "$1/props_cube.sh"; PROPS_TOOLS_DIR="$1"; props_result_file=result.out; has_wfc_pw="$2"; nspin=1; run_cube_props'
                result=subprocess.run(['bash','-e','-c',script,'collector',str(tools),flag],cwd=directory,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertEqual((directory/'result.out').read_text(),'CompareWfc_pass 1\n')

    def test_collector_indented_required_reference_cannot_fall_back(self):
        directory,result=self.run_collector()
        self.assertEqual(result.returncode,0,result.stderr)
        import shutil
        shutil.rmtree(directory/'wfc_pw.ref')
        (directory/'result.ref').write_text('  \tCompareWfc_pass 0\n')
        tools=Path(__file__).resolve().parent
        for flag in ['1','0','']:
            with self.subTest(output_flag=flag):
                (directory/'result.out').unlink()
                script='source "$1/props_cube.sh"; PROPS_TOOLS_DIR="$1"; props_result_file=result.out; has_wfc_pw="$2"; nspin=1; run_cube_props'
                result=subprocess.run(['bash','-e','-c',script,'collector',str(tools),flag],cwd=directory,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertEqual((directory/'result.out').read_text(),'CompareWfc_pass 1\n')

    def test_collector_legacy_cases_keep_their_existing_checks(self):
        directory,result=self.run_collector(False)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual((directory/'result.out').read_text(),'Max_wfc_1 0.5000\nMax_wfc_2 0.5000\nMax_wfc_3 0.5000\n')


if __name__=='__main__':unittest.main()
