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
        # Preserve intentional CRLF fixtures independently of host Git settings.
        self.git("config", "core.autocrlf", "false")
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
        options = {} if "b" in mode else {"encoding": "utf-8", "newline": "\n"}
        with open(target, mode, **options) as handle:
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

    def test_warns_for_default_parameters_added_to_headers(self):
        self.write("source/source_base/defaults.h", "void update_solver(int step = 0);\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_warns_with_success(result, "No new default parameters")

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
        self.assertNotIn("forward declaration", result.stdout.lower())

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

    def test_blocks_new_pragma_once_in_header(self):
        self.write("source/source_base/new_header.h", "#pragma once\n\nint foo();\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "No #pragma once")

    def test_blocks_new_pragma_once_in_cuh(self):
        self.write("source/source_base/kernels/new_header.cuh", "#pragma once\n\nint foo();\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "No #pragma once")

    def test_allows_ifndef_include_guard(self):
        self.write(
            "source/source_base/new_header.h",
            "#ifndef NEW_HEADER_H\n#define NEW_HEADER_H\n\nint foo();\n\n#endif\n",
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("No #pragma once", result.stdout)

    def test_allows_existing_pragma_once_on_unchanged_lines(self):
        self.write("source/source_base/old_header.h", "#pragma once\n\nint foo();\n")
        self.write("source/source_base/CMakeLists.txt", "add_library(old old_header.h)\n")
        self.git("add", ".")
        self.git("commit", "-m", "add header with pragma once")
        base2 = self.git("rev-parse", "HEAD").stdout.strip()
        # Touch only an unrelated line; the #pragma once line is not in the diff.
        self.write("source/source_base/old_header.h", "#pragma once\n\nint foo();\nint bar();\n")
        head = self.commit_change()

        result = self.run_checker("--base", base2, "--head", head)

        self.assertNotIn("No #pragma once", result.stdout)

    def write_exception_event(self, rule, labels, approver="@maintainer"):
        event = self.repo / "event.json"
        body = "\n".join(
            [
                "### Governance Exception",
                f"- Rule: {rule}",
                "- Reason: file is consumed by a generated build list.",
                "- Scope: source/source_base/generated.cpp only.",
                "- User or maintenance risk: low.",
                "- Why the normal rule cannot be followed now: generator lands in a follow-up PR.",
                "- Follow-up cleanup plan: wire the file through CMakeLists.txt in the follow-up.",
                f"- Requested approver: {approver}",
            ]
        )
        event.write_text(
            json.dumps(
                {"pull_request": {"body": body, "labels": [{"name": name} for name in labels]}}
            )
        )
        return event

    def test_approved_exception_downgrades_covered_blocker(self):
        self.write("source/source_base/generated.cpp", "int generated();\n")
        head = self.commit_change()
        event = self.write_exception_event(
            "CMake linkage for new sources", ["governance-exception-approved"]
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_warns_with_success(result, "CMake linkage for new sources")
        self.assertIn("Downgraded from error by the approved governance exception", result.stdout)

    def test_unapproved_exception_keeps_blocker(self):
        self.write("source/source_base/generated.cpp", "int generated();\n")
        head = self.commit_change()
        event = self.write_exception_event("CMake linkage for new sources", [])

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_blocked_by(result, "CMake linkage for new sources")
        self.assertIn("awaiting maintainer approval", result.stdout)

    def test_incomplete_exception_keeps_blocker_even_with_label(self):
        self.write("source/source_base/generated.cpp", "int generated();\n")
        head = self.commit_change()
        event = self.write_exception_event(
            "CMake linkage for new sources", ["governance-exception-approved"], approver=""
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_blocked_by(result, "CMake linkage for new sources")
        self.assertIn("missing fields: Requested approver", result.stdout)

    def test_approved_exception_does_not_downgrade_other_rules(self):
        self.write("source/source_base/generated.cpp", "int generated();\n")
        head = self.commit_change()
        event = self.write_exception_event(
            "Global dependency budget", ["governance-exception-approved"]
        )

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_blocked_by(result, "CMake linkage for new sources")
        self.assertNotIn("Downgraded from error", result.stdout)

    def test_approved_exception_never_downgrades_non_exceptable_rule(self):
        self.write("source/source_base/new_header.h", "#pragma once\n\nint foo();\n")
        head = self.commit_change()
        event = self.write_exception_event("No #pragma once", ["governance-exception-approved"])

        result = self.run_checker("--base", self.base, "--head", head, "--event-path", str(event))

        self.assert_blocked_by(result, "No #pragma once")
        self.assertNotIn("Downgraded from error", result.stdout)

    def test_warns_when_approval_label_has_no_exception_section(self):
        event = self.repo / "event.json"
        event.write_text(
            json.dumps(
                {
                    "pull_request": {
                        "body": "### What's changed?\n- Docs tweak.\n",
                        "labels": [{"name": "governance-exception-approved"}],
                    }
                }
            )
        )

        result = self.run_checker("--event-path", str(event))

        self.assert_warns_with_success(result, "Governance exception record")

    def commit_baseline(self, files):
        for path, content in files.items():
            self.write(path, content)
        self.git("add", ".")
        self.git("commit", "-m", "baseline tests")
        return self.git("rev-parse", "HEAD").stdout.strip()

    def test_allows_conforming_unittests_layout(self):
        self.write("source/source_md/unittests/test_verlet.cpp", "TEST(A, B) {}\n")
        self.write("source/source_md/unittests/test_verlet.cu", "int k();\n")
        self.write("source/source_md/unittests/mock_cell.h", "int mock();\n")
        self.write("source/source_md/unittests/mock_dev.cuh", "int dev();\n")
        self.write("source/source_md/unittests/mock_impl.hpp", "int impl();\n")
        self.write("source/source_md/unittests/CMakeLists.txt", "AddTest(test_verlet)\n")
        self.write("source/source_md/unittests/support/STRU", "ATOMIC_SPECIES\n")
        self.write("source/source_md/unittests/data/test/H2O/run.sh", "echo ok\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("Unit test directory layout", result.stdout)
        self.assertNotIn("Unit test file naming", result.stdout)

    def test_blocks_any_file_added_to_legacy_test_directories(self):
        for path in (
            "source/source_md/test/verlet_test.cpp",
            "source/source_md/test/support/STRU",
            "source/source_md/test/CMakeLists.txt",
            "source/source_io/test_serial/test_reader.cpp",
            "source/source_pw/module_pwdft/tests/helper.h",
            "source/source_md/UnitTests/test_verlet.cpp",
            "source/source_md/test/unittests/test_verlet.cpp",
        ):
            with self.subTest(path=path):
                self.write(path, "int x;\n")
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "Unit test directory layout")
                self.assertIn(path, result.stdout)
                self.git("reset", "--hard", self.base)

    def test_blocks_modified_test_source_in_legacy_directory(self):
        base = self.commit_baseline({"source/source_md/test/verlet_test.cpp": "TEST(A, B) {}\n"})
        self.write("source/source_md/test/verlet_test.cpp", "TEST(A, B) {}\nTEST(A, C) {}\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_blocked_by(result, "Unit test directory layout")
        self.assertIn("Modifies a unit test source", result.stdout)

    def test_allows_cleanup_of_legacy_test_directories(self):
        base = self.commit_baseline(
            {
                "source/source_md/test/verlet_test.cpp": "TEST(A, B) {}\n",
                "source/source_md/test/old_test.cpp": "TEST(O, P) {}\n",
                "source/source_md/test/helper.h": "int helper();\n",
                "source/source_md/test/CMakeLists.txt": "AddTest(verlet_test)\nAddTest(old_test)\n",
            }
        )
        self.git("rm", "-q", "source/source_md/test/old_test.cpp")
        self.write("source/source_md/test/helper.h", "int helper(int);\n")
        self.write("source/source_md/test/CMakeLists.txt", "AddTest(verlet_test)\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertNotIn("Unit test directory layout", result.stdout)

    def test_allows_migrating_legacy_test_into_unittests(self):
        body = "".join(f"TEST(Verlet, Case{i}) {{ EXPECT_EQ({i}, {i}); }}\n" for i in range(20))
        base = self.commit_baseline({"source/source_md/test/verlet_test.cpp": body})
        (self.repo / "source/source_md/unittests").mkdir()
        self.git("mv", "source/source_md/test/verlet_test.cpp", "source/source_md/unittests/test_verlet.cpp")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertNotIn("Unit test directory layout", result.stdout)
        self.assertNotIn("Unit test file naming", result.stdout)

    def test_blocks_code_files_in_unittests_subdirectories(self):
        for path in (
            "source/source_md/unittests/support/helper.h",
            "source/source_md/unittests/gpu/test_verlet.cu",
            "source/source_md/unittests/data/deep/mock.cpp",
            "source/source_md/unittests/support/dev.cuh",
        ):
            with self.subTest(path=path):
                self.write(path, "int x;\n")
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "Unit test directory layout")
                self.assertIn(path, result.stdout)
                self.git("reset", "--hard", self.base)

    def test_blocks_unittests_source_without_test_prefix(self):
        for path in (
            "source/source_md/unittests/verlet_test.cpp",
            "source/source_md/unittests/tmp_mocks.cpp",
            "source/source_md/unittests/kernel.cu",
            "source/source_md/unittests/Test_verlet.cpp",
        ):
            with self.subTest(path=path):
                self.write(path, "int x;\n")
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "Unit test file naming")
                self.git("reset", "--hard", self.base)

    def test_blocks_misnamed_unittests_source_when_modified_or_renamed(self):
        body = "".join(f"int mock_{i}() {{ return {i}; }}\n" for i in range(20))
        base = self.commit_baseline({"source/source_md/unittests/tmp_mocks.cpp": body})
        self.write("source/source_md/unittests/tmp_mocks.cpp", body + "int extra() { return 0; }\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_blocked_by(result, "Unit test file naming")

        self.git("reset", "--hard", base)
        self.git("mv", "source/source_md/unittests/tmp_mocks.cpp", "source/source_md/unittests/mocks.cpp")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_blocked_by(result, "Unit test file naming")

    def test_unittest_layout_ignores_paths_outside_source(self):
        self.write("tests/01_PW/test/run.cpp", "int x;\n")
        self.write("python/pyabacus/tests/helper.cpp", "int y;\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("Unit test directory layout", result.stdout)
        self.assertNotIn("Unit test file naming", result.stdout)

    def test_staged_mode_blocks_new_file_in_legacy_test_directory(self):
        self.write("source/source_md/test/verlet_test.cpp", "TEST(A, B) {}\n")
        self.git("add", ".")

        result = self.run_checker("--staged")

        self.assert_blocked_by(result, "Unit test directory layout")

    def test_blocks_deprecated_cmake_patterns(self):
        cases = (
            ("source/source_base/unittests/CMakeLists.txt", "AddTest(TARGET t LIBS base ${math_libs} SOURCES test_a.cpp)\n"),
            ("source/source_base/unittests/CMakeLists.txt", "remove_definitions(-D__LIBXC)\n"),
            ("source/source_base/CMakeLists.txt", "target_compile_definitions(base PRIVATE USE_LIBXC)\n"),
            ("CMakeLists.txt", "if(ENABLE_RAPIDJSON)\nendif()\n"),
            ("CMakeLists.txt", "if(ENABLE_LIBCOMM)\nendif()\n"),
        )
        for path, content in cases:
            with self.subTest(path=path, content=content.strip()):
                self.write(path, content)
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "CMake deprecated patterns")
                self.git("reset", "--hard", self.base)

    def test_allows_cmake_comments_messages_and_pyabacus_for_deprecated_patterns(self):
        self.write(
            "source/source_base/unittests/CMakeLists.txt",
            "# Replaced ${math_libs}, USE_LIBXC, and remove_definitions(-DUSE_LIBXC).\n"
            "abacus_disable_feature_definitions(__LIBXC)\n"
            "AddTest(TARGET t LIBS base SOURCES test_a.cpp)\n",
        )
        self.write(
            "CMakeLists.txt",
            'message(FATAL_ERROR "ENABLE_RAPIDJSON and ENABLE_LIBCOMM were removed.")\n'
            "message(\n"
            '  WARNING\n'
            '  "USE_LIBXC was renamed to __LIBXC."\n'
            ")\n",
        )
        self.write(
            "python/pyabacus/CMakeLists.txt",
            "target_link_libraries(pyabacus PRIVATE ${math_libs})\n",
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("CMake deprecated patterns", result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_blocks_directory_scoped_and_discovery_cmake_commands_under_source(self):
        for line in (
            "add_definitions(-DFOO)\n",
            "add_compile_definitions(FOO)\n",
            "add_compile_options(-O2)\n",
            "include_directories(${CMAKE_CURRENT_SOURCE_DIR})\n",
            "link_directories(/opt/lib)\n",
            "link_libraries(foo)\n",
            "find_package(Foo REQUIRED)\n",
            "find_library(FOO_LIB foo)\n",
        ):
            with self.subTest(line=line.strip()):
                self.write("source/source_base/CMakeLists.txt", line)
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "CMake directory-scoped commands under source/")
                self.git("reset", "--hard", self.base)

    def test_allows_target_scoped_commands_under_source_and_find_in_root_and_modules(self):
        self.write(
            "source/source_base/CMakeLists.txt",
            "# Do not call include_directories() or find_library() here.\n"
            "add_library(demo OBJECT demo.cpp)\n"
            "target_include_directories(demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})\n"
            "target_compile_definitions(demo PRIVATE DEMO_FLAG)\n"
            "target_compile_options(demo PRIVATE -Wall)\n",
        )
        self.write("CMakeLists.txt", "find_package(Foo REQUIRED)\n")
        self.write("cmake/modules/FindFoo.cmake", "find_library(FOO_LIB NAMES foo)\n")
        self.write("cmake/CollectBuildInfoVars.cmake", "set(ABACUS_FOO_SUPPORT yes)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("CMake directory-scoped commands under source/", result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_blocks_legacy_test_file_staging_in_unittest_cmake(self):
        for content in (
            "install(FILES run.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})\n",
            "file(COPY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})\n",
        ):
            with self.subTest(content=content.strip()):
                self.write("source/source_md/unittests/CMakeLists.txt", content)
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "CMake test file preparation")
                self.git("reset", "--hard", self.base)

    def test_allows_abacus_copy_test_files_and_top_level_source_install(self):
        self.write(
            "source/source_md/unittests/CMakeLists.txt",
            "# Migrated from install() and file(COPY ...)\n"
            "abacus_copy_test_files(support run.sh)\n",
        )
        self.write(
            "source/CMakeLists.txt",
            "install(PROGRAMS ${ABACUS_BIN_PATH} TYPE BIN)\n",
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("CMake test file preparation", result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_warns_when_new_cmake_option_or_module_omits_build_info_update(self):
        self.write("CMakeLists.txt", 'option(ENABLE_FOO "Enable Foo" OFF)\n')
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_warns_with_success(result, "CMake option and dependency sync")

        self.git("reset", "--hard", self.base)
        self.write(
            "CMakeLists.txt",
            'cmake_dependent_option(\n  ENABLE_BAR "Enable Bar" OFF "ENABLE_MPI" OFF)\n',
        )
        self.write("cmake/modules/FindBar.cmake", "find_library(BAR_LIB NAMES bar)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_warns_with_success(result, "CMake option and dependency sync")
        self.assertIn("FindBar.cmake", result.stdout)

    def test_no_option_sync_warning_when_collect_build_info_updated(self):
        self.write("CMakeLists.txt", 'option(ENABLE_FOO "Enable Foo" OFF)\n')
        self.write("cmake/modules/FindFoo.cmake", "find_library(FOO_LIB NAMES foo)\n")
        self.write("cmake/CollectBuildInfoVars.cmake", "set(ABACUS_FOO_SUPPORT yes)\n")
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assertNotIn("CMake option and dependency sync", result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_no_option_sync_warning_when_reordering_existing_options(self):
        base = self.commit_baseline(
            {
                "CMakeLists.txt": (
                    'option(ENABLE_A "A" OFF)\n'
                    'option(ENABLE_B "B" OFF)\n'
                )
            }
        )
        self.write(
            "CMakeLists.txt",
            'option(ENABLE_B "B" OFF)\n'
            'option(ENABLE_A "A" OFF)\n',
        )
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertNotIn("CMake option and dependency sync", result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_staged_mode_blocks_cmake_governance_violations(self):
        self.write("source/source_md/unittests/CMakeLists.txt", "include_directories(/tmp)\n")
        self.git("add", ".")

        result = self.run_checker("--staged")

        self.assert_blocked_by(result, "CMake directory-scoped commands under source/")

    def test_cmake_bracket_comments_preserve_unchanged_openers(self):
        path = "source/source_md/unittests/CMakeLists.txt"
        for equals in ("", "=", "=="):
            with self.subTest(equals=equals):
                opener = "#[" + equals + "[\n"
                closer = "]" + equals + "]\n"
                base = self.commit_baseline({path: opener + "old explanation\n" + closer})
                self.write(
                    path,
                    opener + "include_directories(foo)\n"
                    "file(COPY support DESTINATION build)\n"
                    "target_compile_definitions(foo PRIVATE USE_LIBXC)\n" + closer,
                )
                head = self.commit_change()

                result = self.run_checker("--base", base, "--head", head)

                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertNotIn("CMake", result.stdout)

    def test_cmake_block_comment_close_restores_command_detection(self):
        self.write(
            "source/source_base/CMakeLists.txt",
            "#[=[ ignored\n]=] include_directories(foo)\n",
        )
        head = self.commit_change()

        result = self.run_checker("--base", self.base, "--head", head)

        self.assert_blocked_by(result, "CMake directory-scoped commands under source/")

    def test_cmake_quoted_and_bracket_arguments_are_not_diagnostic_exemptions(self):
        for argument in ('"USE_LIBXC"', '"$<$<BOOL:1>:USE_LIBXC>"', "[=[USE_LIBXC]=]"):
            with self.subTest(argument=argument):
                self.write(
                    "source/source_base/CMakeLists.txt",
                    "target_compile_definitions(base PRIVATE\n  " + argument + "\n)\n",
                )
                head = self.commit_change()

                result = self.run_checker("--base", self.base, "--head", head)

                self.assert_blocked_by(result, "CMake deprecated patterns")
                self.git("reset", "--hard", self.base)

    def test_cmake_message_context_handles_multiline_and_bracket_arguments(self):
        path = "source/source_md/unittests/CMakeLists.txt"
        base = self.commit_baseline({path: 'message(STATUS "old text\nmore text")\n'})
        self.write(
            path,
            'message(STATUS "old text\nUSE_LIBXC ${math_libs}\n'
            'include_directories(foo) # quoted hash\nfile(COPY data)\nmore text")\n'
            'message([==[\nUSE_LIBXC\ninclude_directories(foo)\n]==])\n'
            'set(text "literal \\" # still quoted\ninclude_directories(foo)")\n'
            'set(text [=[\nfile(COPY support)\n]=])\n',
        )
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("CMake", result.stdout)

    def test_cmake_multiline_copy_and_noncopy_control(self):
        path = "source/source_md/unittests/CMakeLists.txt"
        base = self.commit_baseline({path: "file(\n  GLOB files *.txt\n)\n"})
        self.write(path, "file(\n  # stage support files\n  COPY support DESTINATION build\n)\n")
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_blocked_by(result, "CMake test file preparation")
        self.write(path, "file(\n  GLOB files *.cmake\n)\n")
        head = self.commit_change()
        result = self.run_checker("--base", base, "--head", head)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("CMake test file preparation", result.stdout)

    def test_cmake_option_formatting_compares_complete_declaration_names(self):
        multiline = 'cmake_dependent_option(\n  ENABLE_FOO "Foo" OFF "ENABLE_MPI" OFF)\n'
        singleline = 'cmake_dependent_option(ENABLE_FOO "Foo" OFF "ENABLE_MPI" OFF)\n'
        base = self.commit_baseline({"CMakeLists.txt": multiline})
        self.write("CMakeLists.txt", singleline)
        head = self.commit_change()

        for before, after in ((base, head), (head, base)):
            with self.subTest(before=before):
                result = self.run_checker("--base", before, "--head", after)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertNotIn("CMake option and dependency sync", result.stdout)

    def test_cmake_option_rename_on_continuation_line_is_new(self):
        baseline = 'option(\n  ENABLE_OLD "Old" OFF)\n'
        base = self.commit_baseline({"CMakeLists.txt": baseline})
        self.write("CMakeLists.txt", baseline.replace("ENABLE_OLD", "ENABLE_NEW"))
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assert_warns_with_success(result, "CMake option and dependency sync")

    def test_cmake_comment_only_edits_leave_historical_commands_unreported(self):
        path = "source/source_md/unittests/CMakeLists.txt"
        baseline = "include_directories(\n # old note\n foo\n)\nfile(\n # old note\n COPY data DESTINATION build\n)\n"
        base = self.commit_baseline({path: baseline})
        self.write(path, baseline.replace("# old note", "# new note"))
        head = self.commit_change()

        result = self.run_checker("--base", base, "--head", head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("CMake", result.stdout)

    def test_cmake_multiline_arguments_only_report_added_symbol_lines(self):
        path = "source/source_base/CMakeLists.txt"
        for opening, closing in (('[=[', ']=]'), ('"', '"')):
            with self.subTest(opening=opening):
                baseline = "set(NOTES " + opening + "\nUSE_LIBXC is retired\nold wording\n" + closing + ")\n"
                base = self.commit_baseline({path: baseline})
                self.write(path, baseline.replace("old wording", "new wording"))
                head = self.commit_change()
                result = self.run_checker("--base", base, "--head", head)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertNotIn("CMake deprecated patterns", result.stdout)

                self.write(path, baseline.replace("old wording", "ENABLE_LIBCOMM is retired"))
                head = self.commit_change()
                result = self.run_checker("--base", base, "--head", head, "--format", "json")
                findings = json.loads(result.stdout)
                matches = [item for item in findings if item["rule"] == "CMake deprecated patterns"]
                self.assertEqual(len(matches), 1, result.stdout + result.stderr)
                self.assertEqual(matches[0]["line"], 3)
                self.assertIn("ENABLE_LIBCOMM", matches[0]["reason"])

    def test_cmake_staged_context_ignores_unstaged_worktree(self):
        path = "source/source_md/unittests/CMakeLists.txt"
        self.commit_baseline({path: "#[=[\nold note\n]=]\n"})
        self.write(path, "#[=[\ninclude_directories(foo)\n]=]\n")
        self.git("add", ".")
        self.write(path, "include_directories(foo)\n")

        result = self.run_checker("--staged")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("CMake", result.stdout)
        self.write(path, "file(\n COPY data DESTINATION build\n)\n")
        self.git("add", ".")
        self.write(path, "# harmless unstaged replacement\n")
        result = self.run_checker("--staged")
        self.assert_blocked_by(result, "CMake test file preparation")


if __name__ == "__main__":
    unittest.main()
