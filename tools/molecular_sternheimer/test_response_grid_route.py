"""Static production-wiring checks complement the native grid/field regressions."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp'


class MolecularResponseGridRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        text = SOURCE.read_text()
        cls.body = text.split('void run_sternheimer_abacus_chi0_output_impl(',1)[1].split(
            'void run_sternheimer_abacus_chi0_output(',1)[0]

    def test_zero_q_path_reads_response_selection(self):
        for token in ('PARAM.inp.sternheimer_response_nx','PARAM.inp.sternheimer_response_ny',
                      'PARAM.inp.sternheimer_response_nz','make_sternheimer_response_grid(pw_basis,'):
            self.assertTrue(token in self.body,token)
        self.assertLess(self.body.index('run_sternheimer_periodic_lcao_chi0_output('),
                        self.body.index('make_sternheimer_response_grid(pw_basis,'))

    def test_sampling_uses_selected_response_grid(self):
        for token in ('make_sternheimer_fd_full_grid(response_pw_basis)',
                      'make_sternheimer_fd_grid(response_pw_basis)',
                      'build_abfs_ccp_data(ucell, grid_data.grid,',
                      'build_lcao_candidate_grid_functions(ucell, grid_data.grid, lcao_orbitals)'):
            self.assertTrue(token in self.body,token)

    def test_independent_hamiltonian_gathers_correct_spin_and_resamples_projectors(self):
        block = self.body.split('hamiltonians.push_back(',1)[1].split('append_chi0_progress_event(',1)[0]
        for token in ('if (response_grid.independent)',
                      'copy_sternheimer_full_local_potential(potential, pw_basis, response_spin_index)',
                      'restrict_sternheimer_real_field_rectangular(',
                      'restrict_sternheimer_real_field(',
                      'response_grid_data.grid.kpoint = response_grid_kpoint',
                      'make_sternheimer_fd_nonlocal_projector_from_unitcell(',
                      'response_grid_data.grid, response_grid_data.volume_element',
                      'make_sternheimer_fd_hamiltonian_from_local_potential(',
                      'make_sternheimer_fd_full_hamiltonian(',
                      'make_sternheimer_fd_hamiltonian('):
            self.assertTrue(token in block,token)

    def test_report_and_unsupported_route_guard(self):
        for token in ('"sternheimer_response_grid_source "','"pbe_grid "',
                      '"sternheimer_response_grid_requested "','"response_grid "',
                      '(!use_delta_sternheimer || !use_lcao_zero_order || write_siab)'):
            self.assertTrue(token in self.body,token)


if __name__ == '__main__':
    unittest.main()
