"""Static integration guards, not a substitute for native/physical validation."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
CPP = ROOT / 'source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp'


class ParallelRoute(unittest.TestCase):
    def test_parallel_route_owns_mutable_solve_state(self):
        text = CPP.read_text()
        text = text[text.index('void run_sternheimer_weak_q_unit('):]
        text = text[:text.index('void run_sternheimer_weak_response_audit(')]
        for required in ('WEAK_Q_WORKERS', 'WEAK_Q_INNER_THREADS',
                         'run_sternheimer_weak_q_tasks', 'worker_ops.at(worker_index)',
                         'fftw_plan_with_nthreads', 'std::lock_guard<std::mutex>',
                         'parallel_batch', 'worker_tasks', 'benchmark_complete'):
            self.assertIn(required, text)

    def test_reference_import_precedes_response_hash(self):
        text = CPP.read_text()
        text = text[text.index('void run_sternheimer_weak_q_unit('):]
        text = text[:text.index('void run_sternheimer_weak_response_audit(')]
        self.assertIn('WEAK_Q_REFERENCE_IMPORT', text)
        self.assertLess(text.index('WEAK_Q_REFERENCE_IMPORT'), text.index('siab::Sha256 reference_digest'))
        self.assertIn('WEAK_Q_COULOMB_REFERENCE_DIR', text)


if __name__ == '__main__':
    unittest.main()
