"""Static routing checks; numerical spectrum tests run on the cluster."""
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[2] / "source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp"


class SpectrumAuditRoute(unittest.TestCase):
    def test_weak_response_audit_is_opt_in_and_fatal(self):
        text = SOURCE.read_text()
        self.assertIn('env_is_true("ABACUS_STERNHEIMER_WEAK_RESPONSE_AUDIT")', text)
        branch = text[text.rindex("catch (const std::exception& error)"):]
        self.assertIn('env_is_true("ABACUS_STERNHEIMER_WEAK_RESPONSE_AUDIT")', branch)
        self.assertIn("std::exit(EXIT_FAILURE)", branch)

    def test_weak_audit_keeps_fine_vertices_and_complete_target(self):
        text = SOURCE.read_text()
        start = text.index("void run_sternheimer_weak_response_audit(")
        end = text.index("std::vector<std::vector<double>> collect_channel_potentials", start)
        audit = text[start:end]
        self.assertIn("orthonormalize_sternheimer_weak_states_in_place", audit)
        self.assertIn("2.0 * channels[j].potential_r[ir] * source_values[ir]", audit)
        self.assertIn("0.5 * sternheimer_fd_grid_dot(fine_vertices[i]", audit)
        self.assertIn("target_record.unoccupied_eigenvalues.size()", audit)
        self.assertIn("static_cast<std::size_t>(PARAM.globalv.nlocal)", audit)
        self.assertIn("requires nbands=nlocal", audit)
        self.assertIn("blocks->expand_coordinates", audit)
        self.assertNotIn("write_chi0", audit)
        self.assertNotIn("apply_kinetic", audit)

    def test_periodic_excitation_gate_precedes_solver_admission(self):
        text = SOURCE.read_text()
        start = text.index("pair_delta_options.retain_grid_functions = false;")
        stop = text.index('append_chi0_progress_event("delta_subspace_ready"', start)
        branch = text[start:stop]
        self.assertIn("require_positive_periodic_delta_excitations(", branch)
        self.assertIn("source_record.eigenvalues", branch)

    def test_empty_virtual_block_uses_the_same_fatal_gate(self):
        text = SOURCE.read_text()
        start = text.index("pair_delta_options.retain_grid_functions = false;")
        stop = text.index("require_positive_periodic_delta_excitations(", start)
        self.assertNotIn("throw std::runtime_error", text[start:stop])

    def test_excitation_failure_is_fatal_even_for_one_rank(self):
        text = SOURCE.read_text()
        start = text.rindex("catch (const std::exception& error)")
        branch = text[start:]
        self.assertIn("SternheimerExcitationError", branch)
        self.assertIn("std::exit(EXIT_FAILURE)", branch)

    def test_audit_uses_blocked_assembly(self):
        text = SOURCE.read_text()
        start = text.index("void run_sternheimer_spectrum_audit(")
        stop = text.index("std::vector<std::vector<double>> collect_channel_potentials", start)
        self.assertTrue("assemble_delta_sternheimer_grid_matrices_fast(" in text[start:stop])

    def test_audit_precedes_channel_allocation_and_returns(self):
        text = SOURCE.read_text()
        marker = 'if (env_is_true("ABACUS_STERNHEIMER_SPECTRUM_AUDIT")'
        self.assertTrue(marker in text, "Missing opt-in spectrum audit")
        start = text.index(marker)
        stop = text.index('append_chi0_progress_event("abfs_source_ready"', start)
        branch = text[start:stop]
        self.assertIn("run_sternheimer_spectrum_audit", branch)
        self.assertIn("return;", branch)
        self.assertIn("GlobalV::NPROC != 1", branch)
        self.assertIn('|| env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT")', branch)

    def test_opt_in_audit_errors_are_fatal(self):
        text = SOURCE.read_text()
        branch = text[text.rindex("catch (const std::exception& error)"):]
        self.assertIn('env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT")', branch)
        self.assertIn('env_is_true("ABACUS_STERNHEIMER_SPECTRUM_AUDIT")', branch)

    def test_audit_cannot_be_mistaken_for_response_success(self):
        text = SOURCE.read_text()
        for marker in ('"status diagnostic_only', '"physical_result no',
                       '"response_equations 0', '"STERNHEIMER_SPECTRUM_AUDIT.dat"'):
            self.assertTrue(marker in text, "Missing audit marker: " + marker)

    def test_galerkin_audit_uses_fine_potential_and_declares_reference(self):
        text = SOURCE.read_text()
        start = text.index('if (env_is_true("ABACUS_STERNHEIMER_GALERKIN_AUDIT"))')
        stop = text.index("else", start)
        self.assertIn("make_sternheimer_fd_full_grid(pw_basis)", text[start:stop])
        self.assertIn("copy_sternheimer_full_local_potential(potential, pw_basis, 0)", text[start:stop])
        self.assertNotIn("independent_response_full_potential", text[start:stop])
        self.assertIn("reference FullGrid_FD8_not_weak_A", text)
        self.assertIn("full_pc_positivity_certified no", text)

    def test_galerkin_audit_preserves_raw_norms_and_all_states(self):
        text = SOURCE.read_text()
        start = text.index("void run_sternheimer_galerkin_spectrum_audit(")
        stop = text.index("std::vector<std::vector<double>> collect_channel_potentials", start)
        branch = text[start:stop]
        self.assertIn("target_record, 0, true, 0,", branch)
        self.assertIn("raw_norms.at(state)", branch)
        self.assertIn("for (const bool identity : {true, false})", branch)


if __name__ == "__main__":
    unittest.main()
