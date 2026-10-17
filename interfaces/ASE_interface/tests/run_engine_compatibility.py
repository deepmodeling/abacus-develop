#!/usr/bin/env python3
"""Run the installed abacuslite wheel against one explicit ABACUS executable.

The runner copies only this interface harness into a report directory outside
the ABACUS checkout. Each real-engine script then gets its own Python process.
"""

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import zipfile
from email.parser import BytesParser
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CASES = {
    'latest': {'scf.py': 1, 'relax.py': 1, 'md.py': 1, 'native_md.py': 1,
               'band.py': 1, 'magnetic.py': 3, 'socket_handoff.py': 1},
    'legacy': {'scf.py': 1, 'relax.py': 1, 'md.py': 1, 'native_md.py': 1,
               'band.py': 1, 'magnetic.py': 3, 'native_md_generated.py': 1},
}
ASSETS = ('Si_ONCV_PBE-1.0.upf', 'Si_gga_8au_100Ry_2s2p1d.orb',
          'Fe_ONCV_PBE-1.0.upf', 'Fe.upf', 'Fe_gga_6au_100Ry_4s2p2d1f.orb')
PROBE = r'''import json, sys
from importlib import metadata
from pathlib import Path
import abacuslite
import abacuslite.core as core
from abacuslite import Abacus, AbacusProfile
d = metadata.distribution('abacuslite')
p = AbacusProfile(command=sys.argv[1])
c = Abacus(profile=p)
print(json.dumps({'python': sys.executable, 'prefix': sys.prefix,
                  'module': str(Path(abacuslite.__file__).resolve()),
                  'installed': str(Path(d.locate_file('abacuslite/__init__.py')).resolve()),
                  'version': d.version, 'engine_version': p.version(),
                  'reader_selection': ('instance' if hasattr(c.template, '_legacyio')
                                       else 'global'),
                  'legacyio': getattr(c.template, '_legacyio', core.__LEGACYIO__)}))'''


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_cases(fmt, results):
    """Reject omitted, skipped or misreported required scripts."""
    expected = CASES[fmt]
    if [item['script'] for item in results] != list(expected):
        raise ValueError(f'incomplete {fmt} lane: expected {list(expected)}')
    for item in results:
        if (item['returncode'] != 0 or item['timed_out'] or
                item['tests_reported'] != expected[item['script']] or
                item['skipped']):
            raise ValueError(f"required case failed: {item['script']}: {item}")


def capture(command, *, cwd, env, log, timeout):
    start = time.monotonic()
    timed_out = False
    with Path(log).open('w') as stream:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=10)
    return {'returncode': process.returncode, 'timed_out': timed_out,
            'seconds': round(time.monotonic() - start, 3)}


def wheel_runtime(wheel):
    """Read wheel version and package runtime hashes without importing it."""
    with zipfile.ZipFile(wheel) as archive:
        metadata_names = [name for name in archive.namelist()
                          if name.endswith('.dist-info/METADATA')]
        if len(metadata_names) != 1:
            raise ValueError('wheel has no unique METADATA')
        metadata = BytesParser().parsebytes(archive.read(metadata_names[0]))
        if metadata.get('Name', '').lower() != 'abacuslite':
            raise ValueError('wheel is not abacuslite')
        runtime = {name: hashlib.sha256(archive.read(name)).hexdigest()
                   for name in archive.namelist()
                   if name.startswith('abacuslite/') and name.endswith('.py')}
        if 'abacuslite/__init__.py' not in runtime:
            raise ValueError('wheel lacks abacuslite runtime')
    return metadata['Version'], runtime


def harness_revision(repo):
    """Record checkout identity separately from copied script hashes."""
    revision = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=False)
    status = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain'],
                            capture_output=True, text=True, check=False)
    return {'ref': revision.stdout.strip() if revision.returncode == 0 else None,
            'dirty': bool(status.stdout.strip()) if status.returncode == 0 else None}


def attest_install(probe, wheel, repo):
    """Require every wheel runtime module to match installed site-packages."""
    version, runtime = wheel_runtime(wheel)
    imported = Path(probe['module']).resolve()
    installed = Path(probe['installed']).resolve()
    prefix = Path(probe['prefix']).resolve()
    if probe['version'] != version or imported != installed:
        raise ValueError(f'installed wheel identity mismatch: {probe}')
    if not imported.is_relative_to(prefix) or 'site-packages' not in imported.parts:
        raise ValueError(f'not installed in interpreter site-packages: {imported}')
    if imported.is_relative_to(repo):
        raise ValueError(f'ABACUS checkout shadowed installed wheel: {imported}')
    package = imported.parent.parent
    for name, digest in runtime.items():
        path = package / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f'installed runtime differs from supplied wheel: {name}')
    return {'version': version, 'runtime_sha256': runtime}


def required_path(value, *, file=False, directory=False, follow_symlinks=True):
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f'path must be absolute: {value}')
    path = path.resolve() if follow_symlinks else path.absolute()
    if file and not path.is_file():
        raise FileNotFoundError(path)
    if directory and not path.is_dir():
        raise FileNotFoundError(path)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', required=True)
    parser.add_argument('--mpi-command', required=True)
    parser.add_argument('--assets', required=True)
    parser.add_argument('--format', choices=CASES, required=True)
    parser.add_argument('--wheel', required=True)
    parser.add_argument('--report-dir', required=True)
    parser.add_argument('--engine-ref', required=True)
    parser.add_argument('--python', default=sys.executable)
    args = parser.parse_args(argv)
    if not re.fullmatch(r'[0-9a-f]{40}', args.engine_ref):
        parser.error('--engine-ref must be a full lowercase commit SHA')
    engine = required_path(args.engine, file=True)
    # Resolving a venv's python symlink would select its base interpreter.
    python = required_path(args.python, file=True, follow_symlinks=False)
    assets = required_path(args.assets, directory=True)
    wheel = required_path(args.wheel, file=True)
    report = required_path(args.report_dir)
    if report.is_relative_to(REPO) or REPO.is_relative_to(report):
        raise ValueError('report directory must be outside ABACUS checkout')
    if report.exists() and any(report.iterdir()):
        raise ValueError(f'report directory is not empty: {report}')
    missing = [name for name in ASSETS if not (assets / name).is_file()]
    if missing:
        raise FileNotFoundError(f'missing curated assets: {missing}')
    mpi = shlex.split(args.mpi_command)
    if not mpi:
        raise ValueError('empty MPI command')
    report.mkdir(parents=True, exist_ok=True)
    harness = report / 'harness'
    harness.mkdir()
    scripts = set(CASES[args.format]) | {'common.py', 'legacy_facts.py'}
    for name in scripts:
        shutil.copy2(HERE / name, harness / name)
    env = os.environ.copy()
    for key in ('PYTHONPATH', 'PYTHONHOME', 'ASE_ABACUS_SOCKET', 'ASE_ABACUS_RAW_DIR'):
        env.pop(key, None)
    env.update(ABACUS_COMMAND=shlex.join(mpi + [str(engine)]),
               ASE_ABACUS_ASSETS=str(assets), ASE_ABACUS_FORMAT=args.format,
               ASE_ABACUS_KEEP_DIR=str(report / 'failed-runs'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    record = {'schema': 'abacuslite-engine-compatibility-v1', 'format': args.format,
              'engine': str(engine), 'engine_ref_asserted': args.engine_ref,
              'engine_sha256': sha256(engine), 'mpi_command': mpi,
              'harness_checkout': harness_revision(REPO),
              'assets': {name: sha256(assets / name) for name in ASSETS},
              'wheel': str(wheel), 'wheel_sha256': sha256(wheel),
              'python': str(python), 'harness_sha256': sha256(__file__),
              'scripts_sha256': {name: sha256(harness / name) for name in sorted(scripts)},
              'expected_cases': CASES[args.format], 'cases': [], 'all_pass': False}
    output = report / 'evidence.json'

    def save():
        output.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')

    save()
    try:
        probe_result = capture([str(engine), '--version'], cwd=harness, env=env,
                               log=report / 'engine-version.log', timeout=30)
        record['engine_probe'] = probe_result
        record['engine_version_output'] = (report / 'engine-version.log').read_text()
        save()
        if probe_result['returncode'] != 0 or probe_result['timed_out']:
            raise RuntimeError('engine --version failed')
        probe_result = capture([str(python), '-c', PROBE, str(engine)], cwd=harness,
                               env=env, log=report / 'package-probe.log', timeout=45)
        record['package_probe_result'] = probe_result
        if probe_result['returncode'] != 0 or probe_result['timed_out']:
            raise RuntimeError('installed package probe failed')
        probe = json.loads((report / 'package-probe.log').read_text().splitlines()[-1])
        record['package_probe'] = probe
        record['wheel_runtime'] = attest_install(probe, wheel, REPO)
        if probe['legacyio'] is not (args.format == 'legacy'):
            raise ValueError(f'wrong calculator reader for {args.format}: {probe}')
        save()
        for script, expected in CASES[args.format].items():
            name = script.removesuffix('.py')
            raw = report / 'raw' / name
            raw.mkdir(parents=True)
            case_env = dict(env, ASE_ABACUS_RAW_DIR=str(raw))
            log = report / f'{name}.log'
            result = capture([str(python), str(harness / script), '-v'], cwd=harness,
                             env=case_env, log=log, timeout=240)
            content = log.read_text()
            counts = [int(value) for value in re.findall(r'Ran (\d+) tests?', content)]
            result.update(script=script, expected_tests=expected,
                          tests_reported=sum(counts),
                          skipped=bool(re.search(r'\bskipped=\d+\b', content)),
                          log=log.name)
            record['cases'].append(result)
            save()
            print(f'{script}: {result}', flush=True)
        validate_cases(args.format, record['cases'])
        record['all_pass'] = True
    except Exception as exc:
        record['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        record['raw_sha256'] = {str(path.relative_to(report)): sha256(path)
                                for path in sorted((report / 'raw').rglob('*'))
                                if path.is_file()} if (report / 'raw').exists() else {}
        record['dependency_freeze'] = None
        freeze = capture([str(python), '-m', 'pip', 'freeze'], cwd=harness, env=env,
                         log=report / 'pip-freeze.log', timeout=45)
        record['dependency_freeze'] = freeze
        save()
    return 0


if __name__ == '__main__':
    sys.exit(main())
