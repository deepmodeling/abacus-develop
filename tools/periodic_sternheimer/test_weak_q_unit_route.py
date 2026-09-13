"""Static contracts only. Native physics tests must run on the remote host."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
CPP = ROOT / "source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp"
HEADER = ROOT / "source/source_lcao/module_ri/sternheimer_weak_q_unit.h"


class WeakQUnitRoute(unittest.TestCase):
    def route(self):
        text = CPP.read_text()
        self.assertTrue("void run_sternheimer_weak_q_unit(" in text, "q-unit driver is missing")
        start = text.index("void run_sternheimer_weak_q_unit(")
        stop = text.index("void run_sternheimer_weak_response_audit(", start)
        return text[start:stop]

    def test_pure_config_has_no_environment_or_io(self):
        self.assertTrue(HEADER.exists(), "pure q-unit configuration is missing")
        text = HEADER.read_text()
        for forbidden in ("getenv(", "MPI_", "ofstream", "PARAM."):
            self.assertNotIn(forbidden, text)

    def test_only_coulomb_bootstrap_may_generate_frequency_grid(self):
        text = CPP.read_text()
        start = text.index('const bool use_weak_q_unit =')
        guard = text[start:text.index('if (PARAM.inp.nspin != 1)', start)]
        condition = ('(PARAM.inp.sternheimer_frequency_grid_file.empty() '
                     '&& !env_is_true("ABACUS_STERNHEIMER_WEAK_Q_COULOMB_ONLY"))')
        self.assertTrue(condition in ' '.join(guard.split()),
                        'missing fixed grid must reject response units but admit Coulomb-only bootstrap')
        self.assertIn('if (use_weak_q_unit', guard)
        self.assertIn('throw std::invalid_argument', guard)
        route = self.route()
        self.assertIn('"frequency_table " << nfreq', route)
        self.assertIn('for (int f = 0; f < nfreq; ++f)', route)
        self.assertIn('frequency_grid.weights_ha[f]', route)
        self.assertLess(route.index('manifest.close()'), route.index('if (coulomb_only)'))

    def test_opt_in_bypasses_response_symmetry_not_pbe_symmetry(self):
        text = CPP.read_text()
        self.assertTrue('env_is_true("ABACUS_STERNHEIMER_WEAK_Q_UNIT")' in text)
        self.assertTrue('PARAM.inp.symmetry == "1" && !use_weak_q_unit' in text)
        self.assertNotIn('PARAM.inp.symmetry = "-1"', text)
        self.assertIn('env_is_true("ABACUS_STERNHEIMER_WEAK_Q_UNIT")', text[text.rindex("catch (const std::exception& error)"):])

    def test_exact_2d_and_folded_full_bloch_contract(self):
        route = self.route()
        self.assertIn("validate_sternheimer_weak_q_fold", route)
        self.assertIn("pair.reciprocal_shift", route)
        self.assertIn("solve_sternheimer_abf_strict2d_coulomb_in_place", route)
        self.assertNotIn("solve_sternheimer_abf_periodic_full_coulomb", route)
        self.assertNotIn("std::exp", route)
        self.assertIn("fine_grid.grid.kpoint = sternheimer_lcao_grid_kpoint(target_record)", route)

    def test_blocks_outside_bands_and_vertices_outside_frequencies(self):
        route = self.route()
        block = route.index("op->assemble_blocks(")
        band = route.index("for (int ib = unit.band_begin")
        vertex = route.index("projected_vertices.f.resize")
        frequency = route.index("for (int ifrequency = unit.frequency_begin", vertex)
        self.assertLess(block, band)
        self.assertLess(band, vertex)
        self.assertLess(vertex, frequency)
        self.assertIn("target_record, 0, true, 0,", route)
        self.assertNotIn("VBM_only", route)

    def test_weighted_partial_columns_never_reader_or_double_hermitization(self):
        route = self.route()
        for required in ("sternheimer_weak_q_weight", "accumulate_sternheimer_weak_q_column",
                         "for (int sign : {1, -1})", "occupation_weighted yes",
                         "full_q_response no", "physical_result no", "reader_v1 no",
                         "validate_sternheimer_weak_q_residual", "vertex_check",
                         "unit_complete yes", "expected_equations"):
            self.assertIn(required, route)
        for forbidden in ("write_chi0", "symmetrize_chi0", "make_sternheimer_partial_response_record"):
            self.assertNotIn(forbidden, route)

    def test_full_input_manifest_and_identity(self):
        route = self.route()
        for required in ("input_manifest_sha256", "reference_sha256", "auxiliary_sha256",
                         "source_band", "kweight", "occupation", "frequency_table",
                         "weights_ha", "full_kpoints", "reciprocal_shift", "source_k",
                         "band_begin", "band_end", "frequency_begin", "frequency_end"):
            self.assertIn(required, route)

    def test_coulomb_gate_precedes_target_blocks(self):
        route = self.route()
        self.assertIn('ABACUS_STERNHEIMER_WEAK_Q_COULOMB_ONLY', route)
        self.assertLess(route.index('if (coulomb_only)'), route.index('op->assemble_blocks('))
        for required in ('coulomb_integral', 'kernel_matched_to_ewald no', 'sample_density',
                         'coulomb_sample_channels', 'unit.stem(response_plan.iq)'):
            self.assertIn(required, route)

    def test_source_band_order_and_original_energies_are_preserved(self):
        route = self.route()
        source = route[route.index('for (int ib = unit.band_begin'):]
        self.assertIn('source_record.coefficients.at(ib - 1)', source)
        self.assertIn('source_record.eigenvalues.at(ib - 1)', source)
        self.assertNotIn('orthonormalize_sternheimer_weak_states_in_place', source)
        self.assertNotIn('std::sort', source)
        self.assertNotIn('max_element', source)

    def test_cache_opt_in_has_before_after_original_operator_probe(self):
        route = self.route()
        self.assertIn('ABACUS_STERNHEIMER_WEAK_EXACT_CACHE', route)
        self.assertLess(route.index('op->apply(probe, uncached)'), route.index('enable_exact_apply_cache'))
        self.assertLess(route.index('enable_exact_apply_cache'), route.index('op->apply(probe, cached)'))
        self.assertIn('cache_error > 1e-10', route)
        self.assertIn('cache_probe_is_not_operator_proof yes', route)

    def test_density_support_and_same_q_ewald_identity_are_not_silently_assumed(self):
        route = self.route()
        self.assertIn('sternheimer_weak_q_radial_support', route)
        self.assertIn('sternheimer_weak_q_z_support_contained', route)
        self.assertIn('coulomb_reference_file', route)
        self.assertIn('find_coulomb_v1_rank_files(response_plan.iq', route)
        self.assertLess(route.index('if (!z_support_contained)'), route.index('op->assemble_blocks('))

    def test_fine_checks_only_first_last_owned_column_all_block_columns_remain(self):
        route = self.route()
        self.assertTrue('if (owned == 0 || owned + 1 == owned_columns.size())' in route)
        selected = route.index('if (owned == 0 || owned + 1 == owned_columns.size())')
        self.assertLess(selected, route.index('worker_ops.at(worker_index)->lift(expansion.e, fine_response)'))
        self.assertLess(selected, route.index('const Vector diagonal_vertex'))
        self.assertLess(route.index("BlasConnector::gemv('C'"), selected)
        self.assertIn('fine_checked_equations', route)
        self.assertIn('fine_vertex_checks', route)


if __name__ == "__main__":
    unittest.main()
