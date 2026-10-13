"""Guard the optional libvdwxc integration contract."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


class LibvdwxcContract(unittest.TestCase):
    def test_optional_backend_has_a_separate_build_switch_and_adapter(self):
        cmake = (ROOT / "CMakeLists.txt").read_text()
        module_cmake = (ROOT / "source/source_hamilt/module_vdw/CMakeLists.txt").read_text()
        adapter = ROOT / "source/source_hamilt/module_vdw/nonlocal_xc_libvdwxc.cpp"
        self.assertIn("ENABLE_LIBVDWXC", cmake)
        self.assertIn("find_package(libvdwxc CONFIG REQUIRED)", cmake)
        self.assertIn("libvdwxc::vdwxc", cmake)
        self.assertIn("ENABLE_LIBVDWXC", module_cmake)
        self.assertIn("nonlocal_xc_libvdwxc.cpp", module_cmake)
        self.assertTrue(adapter.is_file())
        text = adapter.read_text()
        self.assertIn("vdwxc_calculate", text)
        self.assertIn("vdwxc_set_unit_cell", text)
        self.assertIn("vdwxc_init_serial", text)
        self.assertIn("vdwxc_finalize", text)
        adapter_header = (ROOT / "source/source_hamilt/module_vdw/nonlocal_xc_libvdwxc.h").read_text()
        self.assertLess(adapter_header.index("#include <complex.h>"), adapter_header.index("#include <vdwxc.h>"))

    def test_first_pr_does_not_bundle_a_second_rvv10_kernel(self):
        module_cmake = (ROOT / "source/source_hamilt/module_vdw/CMakeLists.txt").read_text()
        self.assertNotIn("nonlocal_xc_rvv10.cpp", module_cmake)
        self.assertNotIn("nonlocal_xc_rvv10_pw.cpp", module_cmake)
        for name in (
            "nonlocal_xc_rvv10.cpp",
            "nonlocal_xc_rvv10.h",
            "nonlocal_xc_rvv10_pw.cpp",
            "nonlocal_xc_rvv10_pw.h",
        ):
            self.assertFalse((ROOT / "source/source_hamilt/module_vdw" / name).exists(), name)

    def test_configuration_checks_the_complete_serial_rvv10_api(self):
        cmake = (ROOT / "CMakeLists.txt").read_text()
        for symbol in ("VDWXC_RVV10", "vdwxc_new", "vdwxc_set_rvv10_parameters", "vdwxc_set_unit_cell",
                       "vdwxc_init_serial", "vdwxc_calculate", "vdwxc_finalize"):
            self.assertIn(symbol, cmake)
        self.assertIn("NO_DEFAULT_PATH", cmake)

    def test_input_path_requires_both_semilocal_and_nonlocal_dependencies(self):
        source = (ROOT / "source/source_io/module_parameter/read_inp_estruc.cpp").read_text()
        cmake = (ROOT / "CMakeLists.txt").read_text()
        self.assertIn("#if !defined(__LIBXC) || !defined(__LIBVDWXC)", source)
        self.assertIn("ENABLE_LIBVDWXC currently requires ENABLE_LIBXC=ON", cmake)

    def test_docs_do_not_describe_a_bundled_native_kernel(self):
        docs = (ROOT / "docs/developers_guide/nonlocal_xc.md").read_text()
        self.assertNotIn("existing native evaluator remains available", docs)
        self.assertIn("does not carry a second rVV10 kernel", docs)

    def test_build_info_reports_the_optional_backend(self):
        build_info = (ROOT / "source/source_io/build_info.h.in").read_text()
        parse_args = (ROOT / "source/source_io/parse_args.cpp").read_text()
        collect_build_info = (ROOT / "cmake/CollectBuildInfoVars.cmake").read_text()
        make_build_info = (ROOT / "generate_build_info.sh").read_text()
        self.assertIn("ABACUS_LIBVDWXC_VERSION", build_info)
        self.assertIn('print_info("libvdwxc Support", ABACUS_LIBVDWXC_VERSION)', parse_args)
        self.assertIn('libvdwxc       = ${ABACUS_LIBVDWXC_VERSION}', collect_build_info)
        self.assertIn('@ABACUS_LIBVDWXC_VERSION@', make_build_info)


if __name__ == "__main__":
    unittest.main()
