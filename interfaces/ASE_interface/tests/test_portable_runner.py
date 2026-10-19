"""Offline contract checks for the portable real-engine lane."""

import os
import signal
import sys
import tempfile
import time
import unittest
import zipfile
from pathlib import Path

from run_engine_compatibility import (CASES, attest_install, capture,
                                      required_path, validate_cases)


class PortableRunnerTest(unittest.TestCase):
    def test_venv_python_symlink_is_not_resolved_to_base_interpreter(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / 'base-python'
            binary.touch()
            venv_python = root / 'venv-python'
            venv_python.symlink_to(binary)
            self.assertEqual(required_path(str(venv_python), file=True,
                                           follow_symlinks=False), venv_python)

    def test_required_case_set_rejects_missing_failure_skip_and_wrong_count(self):
        cases = [{'script': name, 'returncode': 0, 'timed_out': False,
                  'tests_reported': count, 'skipped': False}
                 for name, count in CASES['legacy'].items()]
        validate_cases('legacy', cases)
        for changed in (cases[:-1],):
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                validate_cases('legacy', changed)
        for change in ({'returncode': 1}, {'timed_out': True},
                       {'tests_reported': 0}, {'skipped': True}):
            altered = [dict(item) for item in cases]
            altered[0].update(change)
            with self.assertRaisesRegex(ValueError, 'required case failed'):
                validate_cases('legacy', altered)

    def test_installed_runtime_bytes_must_match_supplied_wheel(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            site = root / 'venv' / 'lib' / 'python3.10' / 'site-packages'
            package = site / 'abacuslite'
            package.mkdir(parents=True)
            imported = package / '__init__.py'
            imported.write_text('original')
            wheel = root / 'abacuslite-1.1.0-py3-none-any.whl'
            with zipfile.ZipFile(wheel, 'w') as archive:
                archive.writestr('abacuslite/__init__.py', 'original')
                archive.writestr('abacuslite-1.1.0.dist-info/METADATA',
                                 'Metadata-Version: 2.1\nName: abacuslite\nVersion: 1.1.0\n')
            probe = {'module': str(imported), 'installed': str(imported),
                     'prefix': str(root / 'venv'), 'version': '1.1.0'}
            attest_install(probe, wheel, root / 'checkout')
            imported.write_text('mutated')
            with self.assertRaisesRegex(ValueError, 'runtime differs'):
                attest_install(probe, wheel, root / 'checkout')
            imported.write_text('original')
            with self.assertRaisesRegex(ValueError, 'shadowed'):
                attest_install(probe, wheel, root)

    def test_hung_case_has_bounded_process_group_timeout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = capture(['python3', '-c', 'import time; time.sleep(60)'],
                             cwd=root, env=None, log=root / 'hung.log', timeout=0.1)
            self.assertTrue(result['timed_out'])
            self.assertNotEqual(result['returncode'], 0)
            self.assertLess(result['seconds'], 3)

    def test_timeout_kills_sigterm_ignoring_child_after_leader_exits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pid_file = root / 'child.pid'
            child_code = ('import os,signal,sys,time; '
                          'signal.signal(signal.SIGTERM, signal.SIG_IGN); '
                          'open(sys.argv[1], "w").write(str(os.getpid())); '
                          'time.sleep(60)')
            leader_code = ('import subprocess,sys,time; '
                           'subprocess.Popen([sys.executable, "-c", '
                           'sys.argv[1], sys.argv[2]]); time.sleep(60)')
            child_pid = None
            try:
                result = capture([sys.executable, '-c', leader_code,
                                  child_code, str(pid_file)],
                                 cwd=root, env=None, log=root / 'hung-group.log',
                                 timeout=1)
                self.assertTrue(pid_file.is_file(), 'child probe did not start')
                child_pid = int(pid_file.read_text())
                self.assertTrue(result['timed_out'])
                self.assertLess(result['seconds'], 4)
                for _ in range(50):
                    try:
                        state = Path(f'/proc/{child_pid}/stat').read_text().split()[2]
                    except FileNotFoundError:
                        break
                    if state == 'Z':
                        break
                    time.sleep(0.02)
                else:
                    self.fail(f'timeout left SIGTERM-ignoring child {child_pid} alive')
            finally:
                if child_pid is None and pid_file.is_file():
                    child_pid = int(pid_file.read_text())
                if child_pid is not None:
                    try:
                        child_group = os.getpgid(child_pid)
                        if child_group != os.getpgrp():
                            os.killpg(child_group, signal.SIGKILL)
                        else:
                            os.kill(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

    def test_nonzero_leader_exit_reclaims_remaining_child(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pid_file = root / 'child.pid'
            child_code = ('import os,signal,sys,time; '
                          'signal.signal(signal.SIGTERM, signal.SIG_IGN); '
                          'open(sys.argv[1], "w").write(str(os.getpid())); '
                          'time.sleep(60)')
            leader_code = ('import os,sys,time,subprocess; '
                           'subprocess.Popen([sys.executable, "-c", '
                           'sys.argv[1], sys.argv[2]]); '
                           'deadline=time.time()+10\n'
                           'while not os.path.exists(sys.argv[2]) and '
                           'time.time()<deadline:\n'
                           '    time.sleep(0.01)\n'
                           'sys.exit(3)')
            child_pid = None
            try:
                result = capture([sys.executable, '-c', leader_code,
                                  child_code, str(pid_file)],
                                 cwd=root, env=None,
                                 log=root / 'nonzero-leader.log', timeout=30)
                self.assertFalse(result['timed_out'])
                self.assertEqual(result['returncode'], 3)
                self.assertTrue(pid_file.is_file(), 'child probe did not start')
                child_pid = int(pid_file.read_text())
                self.assertLess(result['seconds'], 5)
                for _ in range(50):
                    try:
                        state = Path(f'/proc/{child_pid}/stat').read_text().split()[2]
                    except FileNotFoundError:
                        break
                    if state == 'Z':
                        break
                    time.sleep(0.02)
                else:
                    self.fail(f'nonzero leader exit left child {child_pid} alive')
            finally:
                if child_pid is None and pid_file.is_file():
                    child_pid = int(pid_file.read_text())
                if child_pid is not None:
                    try:
                        child_group = os.getpgid(child_pid)
                        if child_group != os.getpgrp():
                            os.killpg(child_group, signal.SIGKILL)
                        else:
                            os.kill(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass


if __name__ == '__main__':
    unittest.main()
