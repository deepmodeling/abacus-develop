import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKER = REPO_ROOT / "tools" / "03_code_analysis" / "agent_governance_check.py"


class AgentGovernanceCheckTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name)
        self.git("init")
        self.git("config", "user.email", "agent-governance@example.com")
        self.git("config", "user.name", "Agent Governance Test")
        self.write("README.md", "baseline\n")
        self.write(
            "source/source_io/module_parameter/read_input_item_model.cpp",
            'Input_Item item("old_switch");\n',
        )
        self.git("add", ".")
        self.git("commit", "-m", "baseline")
        self.base = self.git("rev-parse", "HEAD").stdout.strip()

    def tearDown(self):
        self.tmp.cleanup()

    def git(self, *args):
        return subprocess.run(
            ["git", *args],
            cwd=self.repo,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )

    def write(self, path, content, mode="w"):
        target = self.repo / path
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, mode) as handle:
            handle.write(content)

    def commit_change(self):
        self.git("add", ".")
        self.git("commit", "-m", "change")
        return self.git("rev-parse", "HEAD").stdout.strip()

    def run_checker(self, *args):
        return subprocess.run(
            [sys.executable, str(CHECKER), *args],
            cwd=self.repo,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def assert_blocked_by(self, result, rule):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(rule, result.stdout)

    def assert_warns_with_success(self, result, rule):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(rule, result.stdout)

    def test_detects_crlf_in_changed_text_file(self):
        self.write("source/source_base/crlf.cpp", b"int x = 1;\r\n", mode="wb")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "LF line endings")

    def test_allows_escaped_crlf_text_in_changed_text_file(self):
        self.write("source/source_base/escaped.cpp", 'const char* eol = "\\r\\n";\n')
        self.write("source/source_base/CMakeLists.txt", "add_library(escaped escaped.cpp)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_allows_crlf_in_windows_scripts(self):
        self.write("tools/install.bat", b"echo ok\r\n", mode="wb")
        self.write("tools/install.cmd", b"echo ok\r\n", mode="wb")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_blocks_when_global_dependency_budget_increases(self):
        self.write("source/source_base/global.cpp", "int n = GlobalV::NPROC + PARAM.inp.nbands;\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(global global.cpp)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "Global dependency budget")
        self.assertIn("net_delta=2", result.stdout)

    def test_warns_when_global_dependency_usage_is_rebalanced(self):
        self.write("source/source_base/global.cpp", "int old_n = PARAM.inp.nbands;\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(global global.cpp)\n")
        self.git("add", ".")
        self.git("commit", "-m", "add baseline global usage")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write("source/source_base/global.cpp", "int moved_n = GlobalV::NPROC;\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_warns_with_success(result, "Global dependency budget")
        self.assertIn("net_delta=0", result.stdout)

    def test_allows_global_dependency_budget_reduction(self):
        self.write("source/source_base/global.cpp", "int old_n = PARAM.inp.nbands;\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(global global.cpp)\n")
        self.git("add", ".")
        self.git("commit", "-m", "add baseline global usage")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write("source/source_base/global.cpp", "int old_n = 0;\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("Global dependency budget", result.stdout)

    def test_allows_global_names_in_documentation(self):
        self.write("docs/governance-notes.md", "Mention GlobalV::NPROC and PARAM.inp in documentation.\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_allows_global_names_in_governance_checker_tests(self):
        self.write("tools/03_code_analysis/checker_fixture.py", 'pattern = "GlobalV::NPROC and PARAM.inp"\n')
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_blocks_default_parameters_added_to_headers(self):
        self.write("source/source_base/defaults.h", "void update_solver(int step = 0);\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "No new default parameters")

    def test_ignores_crlf_to_lf_only_changes_for_semantic_added_lines(self):
        self.write("source/source_base/defaults.h", b"void update_solver(int step = 0);\r\n", mode="wb")
        self.git("add", ".")
        self.git("commit", "-m", "add crlf header")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write("source/source_base/defaults.h", "void update_solver(int step = 0);\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("No new default parameters", result.stdout)

    def test_staged_mode_ignores_crlf_to_lf_only_semantic_added_lines(self):
        self.write("source/source_base/defaults.h", b"void update_solver(int step = 0);\r\n", mode="wb")
        self.git("add", ".")
        self.git("commit", "-m", "add crlf header")
        self.write("source/source_base/defaults.h", "void update_solver(int step = 0);\n")
        self.git("add", ".")

        result = self.run_checker("--staged")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("No new default parameters", result.stdout)

    def test_allows_for_loop_initializer_in_header(self):
        self.write(
            "source/source_base/loop_header.h",
            "inline int sum(int n) {\n"
            "    int total = 0;\n"
            "    for (int i = 0; i < n; ++i) {\n"
            "        total += i;\n"
            "    }\n"
            "    return total;\n"
            "}\n",
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("No new default parameters", result.stdout)

    def test_warns_but_does_not_block_for_new_hpp_files(self):
        self.write("source/source_base/detail.hpp", "inline int value() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "# listed elsewhere\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Avoid new .hpp propagation", result.stdout)

    def test_blocks_new_source_file_without_cmake_linkage(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "CMake linkage for new sources")

    def test_blocks_input_parameter_changes_without_docs_linkage(self):
        self.write(
            "source/source_io/module_parameter/read_input_item_model.cpp",
            'Input_Item item("new_switch");\nitem.default_value = "0";\n',
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "INPUT parameter documentation linkage")

    def test_allows_parameter_file_comment_only_change_without_docs(self):
        self.write(
            "source/source_io/module_parameter/read_input_item_model.cpp",
            'Input_Item item("old_switch");\n// Keep legacy input switch documented nearby.\n',
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_allows_input_item_test_fixture_without_docs(self):
        self.write("tools/03_code_analysis/input_fixture.py", 'fixture = "Input_Item item(\\"old_switch\\");"\n')
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_allows_input_parameter_changes_with_required_docs(self):
        self.write(
            "source/source_io/module_parameter/read_input_item_model.cpp",
            'Input_Item item("new_switch");\nitem.default_value = "0";\n',
        )
        self.write("docs/parameters.yaml", "parameters: []\n")
        self.write("docs/advanced/input_files/input-main.md", "# INPUT\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_blocks_input_parameter_change_when_required_docs_are_deleted(self):
        self.write("docs/parameters.yaml", "parameters: []\n")
        self.write("docs/advanced/input_files/input-main.md", "# INPUT\n")
        self.git("add", ".")
        self.git("commit", "-m", "add input docs")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write(
            "source/source_io/module_parameter/read_input_item_model.cpp",
            'Input_Item item("new_switch");\nitem.default_value = "0";\n',
        )
        (self.repo / "docs" / "parameters.yaml").unlink()
        (self.repo / "docs" / "advanced" / "input_files" / "input-main.md").unlink()
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_blocked_by(result, "INPUT parameter documentation linkage")

    def test_warns_for_unfilled_pr_template_fields_from_event_payload(self):
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nFix #...\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "- A unit test is added for each new feature or bug fix.\n"
                    }
                }
            )
        )

        result = self.run_checker("--event-path", str(event))

        self.assert_warns_with_success(result, "PR metadata completeness")

    def test_warns_for_empty_pr_template_from_event_payload(self):
        for body in ("", None):
            with self.subTest(body=body):
                event = self.repo / "event.json"
                event.write_text(json.dumps({"pull_request": {"body": body}}))

                result = self.run_checker("--event-path", str(event))

                self.assert_warns_with_success(result, "PR metadata completeness")

    def test_warns_for_missing_pr_body_from_event_payload(self):
        event = self.repo / "event.json"
        event.write_text(json.dumps({"pull_request": {}}))

        result = self.run_checker("--event-path", str(event))

        self.assert_warns_with_success(result, "PR metadata completeness")

    def test_skips_pr_metadata_when_event_payload_is_not_a_pull_request(self):
        event = self.repo / "event.json"
        event.write_text(json.dumps({"workflow_run": {"name": "Agent Governance"}}))

        result = self.run_checker("--event-path", str(event))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("PR metadata completeness", result.stdout)

    def test_accepts_core_pr_template_fields_from_event_payload(self):
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue; governance bootstrap.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Ran python3 -m unittest tools/03_code_analysis/test_agent_governance_check.py.\n\n"
                        "### What's changed?\n"
                        "Adds governance checks only; no runtime behavior change.\n"
                    }
                }
            )
        )

        result = self.run_checker("--event-path", str(event))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("PR metadata completeness", result.stdout)

    def test_accepts_reminder_style_pr_template_from_event_payload(self):
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Reminder\n"
                        "- [ ] I have read `AGENTS.md` and `docs/developers_guide/agent_governance.md`.\n"
                        "- [ ] I have linked an issue or explained why this PR does not need one.\n"
                        "- [ ] I have added adequate unit tests and/or case tests, or explained why not.\n"
                        "- [ ] I have listed the exact verification commands run and their results.\n"
                        "- [ ] I have described user-visible behavior changes, including INPUT parameter changes.\n"
                        "- [ ] I have explained core-module impact for ESolver, HSolver, ElecState, Hamilt, Operator, Psi, or other `source/` changes.\n"
                        "- [ ] I have requested any needed governance exception below.\n\n"
                        "### Linked Issue\nNo issue; governance bootstrap.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Ran python3 -m unittest tools/03_code_analysis/test_agent_governance_check.py.\n\n"
                        "### What's changed?\n"
                        "Adds governance checks only; no runtime behavior change.\n\n"
                        "### Governance Notes\n"
                        "No INPUT, core module, or exception notes.\n"
                    }
                }
            )
        )

        result = self.run_checker("--event-path", str(event))

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("PR metadata completeness", result.stdout)

    def test_warns_for_source_change_without_test_evidence(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(new_feature new_feature.cpp)\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Not filled yet.\n\n"
                        "### What's changed?\n"
                        "Adds a source file.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base helper only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "Test evidence review")

    def test_warns_when_pr_says_no_tests_were_run(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(new_feature new_feature.cpp)\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "No tests were run.\n\n"
                        "### What's changed?\n"
                        "Adds a source file.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base helper only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "Test evidence review")

    def test_warns_when_pr_says_no_tests_added_yet(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(new_feature new_feature.cpp)\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "No tests added yet.\n\n"
                        "### What's changed?\n"
                        "Adds a source file.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base helper only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "Test evidence review")

    def test_accepts_explicit_no_test_rationale(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(new_feature new_feature.cpp)\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Tests not required because documentation only.\n\n"
                        "### What's changed?\n"
                        "Adds a source file.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base helper only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assertNotIn("Test evidence review", result.stdout)

    def test_warns_for_source_change_without_docs_or_no_docs_reason(self):
        self.write("source/source_base/new_feature.cpp", "int new_feature() { return 1; }\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(new_feature new_feature.cpp)\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Ran focused unit tests.\n\n"
                        "### What's changed?\n"
                        "Adds a source file.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base helper only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "Documentation sync review")

    def test_warns_for_source_header_change_without_test_evidence(self):
        self.write("source/source_base/api.h", "void api();\n")
        self.git("add", ".")
        self.git("commit", "-m", "add header")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write("source/source_base/api.h", "#include <vector>\nvoid api();\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_warns_with_success(result, "Test evidence review")

    def test_warns_for_source_header_change_without_docs_or_no_docs_reason(self):
        self.write("source/source_base/api.h", "void api();\n")
        self.git("add", ".")
        self.git("commit", "-m", "add header")
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.write("source/source_base/api.h", "#include <vector>\nvoid api();\n")
        head = self.commit_change()
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### Linked Issue\nNo issue.\n\n"
                        "### Unit Tests and/or Case Tests for my changes\n"
                        "Ran focused unit tests.\n\n"
                        "### What's changed?\n"
                        "Updates a source header.\n\n"
                        "### Governance Checklist\n"
                        "Reviewed.\n\n"
                        "### INPUT Parameter Changes\n"
                        "No INPUT parameter changes.\n\n"
                        "### Core Module Impact\n"
                        "source/source_base API only.\n\n"
                        "### Governance Exception\n"
                        "No exceptions requested.\n"
                    }
                }
            )
        )

        result = self.run_checker("--base", base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "Documentation sync review")

    def test_warns_for_new_header_include(self):
        self.write("source/source_base/include_growth.h", "#include <vector>\nclass IncludeGrowth {};\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_warns_with_success(result, "Header dependency review")

    def test_merge_base_scoped_comparison_excludes_base_branch_only_changes(self):
        self.write("source/source_base/api.h", "void update_solver(int step = 0);\n")
        self.git("add", ".")
        self.git("commit", "-m", "add legacy default")
        base_branch = self.git("branch", "--show-current").stdout.strip()
        merge_base = self.git("rev-parse", "HEAD").stdout.strip()
        self.git("checkout", "-b", "feature")
        self.write("docs/feature.md", "feature docs\n")
        head = self.commit_change()
        self.git("checkout", base_branch)
        self.write("source/source_base/api.h", "void update_solver(int step);\n")
        base_tip = self.commit_change()

        base_tip_result = self.run_checker("--base", base_tip, "--head", head)
        merge_base_result = self.run_checker("--base", merge_base, "--head", head)

        self.assertIn("No new default parameters", base_tip_result.stdout)
        self.assertNotIn("No new default parameters", merge_base_result.stdout)
        self.assertEqual(merge_base_result.returncode, 0, merge_base_result.stdout + merge_base_result.stderr)

    def test_blocks_new_heterogeneous_file_without_cmake_linkage(self):
        self.write("source/module_hamilt/kernels/new_kernel.cu", "__global__ void k() {}\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "CMake linkage for heterogeneous sources")

    def test_blocks_new_cuh_file_without_cmake_linkage(self):
        self.write("source/module_hamilt/kernels/new_kernel.cuh", "__device__ int k();\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "CMake linkage for heterogeneous sources")

    def test_warns_for_heterogeneous_file_without_test_evidence(self):
        self.write("source/module_hamilt/kernels/new_kernel.cu", "__global__ void k() {}\n")
        self.write("source/module_hamilt/CMakeLists.txt", "add_library(new_kernel kernels/new_kernel.cu)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_warns_with_success(result, "Heterogeneous test evidence review")

    def test_staged_mode_blocks_global_dependency_budget_increase(self):
        self.write("source/source_base/staged.cpp", "int n = GlobalC::ucell.nat;\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(staged staged.cpp)\n")
        self.git("add", ".")

        result = self.run_checker("--staged")

        self.assert_blocked_by(result, "Global dependency budget")

    def test_rejects_staged_with_base_head(self):
        result = self.run_checker("--staged", "--base", self.base, "--head", self.base)

        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("--staged cannot be combined with --base/--head", result.stderr)


if __name__ == "__main__":
    unittest.main()
