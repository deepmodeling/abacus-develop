"""Source contract only; native/numerical verification belongs on the remote host."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp"


class WeakMatrixAuditRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = SOURCE.read_text()
        start = cls.text.index("void run_sternheimer_weak_response_audit(")
        stop = cls.text.index("std::vector<std::vector<double>> collect_channel_potentials", start)
        cls.audit = cls.text[start:stop]

    def test_new_switch_enters_existing_route_and_fails_fatally(self):
        switch = 'env_is_true("ABACUS_STERNHEIMER_WEAK_MATRIX_AUDIT")'
        self.assertIn(switch, self.audit)
        route = self.text[self.text.index('if (env_is_true("ABACUS_STERNHEIMER_WEAK_RESPONSE_AUDIT")'):
                          self.text.index('if (env_is_true("ABACUS_STERNHEIMER_SPECTRUM_AUDIT")')]
        self.assertIn(switch, route)
        for token in ("GlobalV::NPROC != 1", "response_plan.kq_pairs.at(0)",
                      "frequency_grid.omega_ha.at(0)", "return;"):
            self.assertIn(token, route)
        self.assertIn(switch, self.text[self.text.rindex("catch (const std::exception& error)"):])

    def test_partial_output_contract(self):
        for token in ("matrix_audit yes", "shard_index ", "shard_count ", "auxiliary_channels ",
                      "auxiliary_channel ", "qpoint ", "response_matrix partial",
                      "global_matrix_checks deferred_external_merge", "source_band_selection VBM_only",
                      "status diagnostic_only", "full_q_response no", "physical_result no",
                      "WEAK_AUDIT_SHARD_INDEX", "WEAK_AUDIT_SHARD_COUNT"):
            self.assertIn(token, self.audit)
        for field in ("atom_index", "atom_local_index", "type_index", "angular_momentum",
                      "magnetic_index", "radial_index", "label"):
            self.assertIn(field, self.audit)
        self.assertNotIn("write_chi0", self.audit)
        self.assertNotIn("MPI_", self.audit)

    def test_matrix_mode_streams_then_caches_projected_vertices_only(self):
        self.assertIn("describe_sternheimer_abf_grid_channels", self.audit)
        self.assertIn("sample_fine_vertex", self.audit)
        self.assertIn("SternheimerWeakAuditShard::parse", self.audit)
        self.assertIn("owned_columns", self.audit)
        self.assertIn("projected_vertices", self.audit)
        self.assertIn("fine_grid.grid, qpoint, 1,", self.audit)
        self.assertNotIn("std::vector<Vector> fine_vertices", self.audit)
        self.assertNotIn("std::vector<Vector> fine_potentials", self.audit)

    def test_projected_contraction_all_rows_fine_check_selected_rows(self):
        self.assertIn("i == j || i == 0", self.audit)
        self.assertIn("sample_fine_vertex(j)", self.audit)
        self.assertIn("sample_fine_vertex(0)", self.audit)
        self.assertIn("block_value", self.audit)
        self.assertIn("vertex_check ", self.audit)
        self.assertIn("response_equations ", self.audit)
        self.assertIn("2 * owned_columns.size()", self.audit)

    def test_original_four_equation_and_full_matrix_gates_remain(self):
        for token in ("std::array<Vector, 2> fine_vertices", "std::array<Complex, 4> response",
                      "response_hermitian_error", "response_max_eigenvalue",
                      "hermitian_error > 1e-7", "maximum_eigenvalue > 1e-8",
                      "auxiliary_channel_selection first_two_only"):
            self.assertIn(token, self.audit)

    def test_complete_input_and_nonfolded_guards_remain(self):
        for token in ("requires nbands=nlocal", "target_record.unoccupied_coefficients.size()",
                      "does not yet admit BZ-folded k pairs", "target_record, 0, true, 0,",
                      "source_record.coefficients.size()", "fine_potential.size()"):
            self.assertIn(token, self.audit)


if __name__ == "__main__":
    unittest.main()
