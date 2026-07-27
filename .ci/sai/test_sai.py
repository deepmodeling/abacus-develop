import argparse
import io
import importlib.util
import json
import subprocess
import sys
import tempfile
import tarfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
SPEC = importlib.util.spec_from_file_location("sai", ROOT / "sai.py")
assert SPEC and SPEC.loader
sai = importlib.util.module_from_spec(SPEC)
sys.modules["sai"] = sai
SPEC.loader.exec_module(sai)
import slurm  # noqa: E402


def valid_result():
    config = sai.load_config()
    components = [sai._component("build", "Compile", "PASS")]
    components.extend(sai._component(name, profile.label, "PASS") for name, profile in config.resources.items())
    rows = [
        sai._result_row(case, "PASS", exit_code=0)
        for resource in config.resources
        for case in config.cases if case.resource == resource
    ]
    return {
        "protocol": 1, "total": len(rows), "passed": len(rows),
        "failed": 0, "infrastructure": 0, "components": components,
        "cases": rows,
    }


class ConfigTests(unittest.TestCase):
    def test_current_matrix_is_loaded_from_ini(self):
        config = sai.load_config()
        self.assertEqual(len(config.cases), 49)
        self.assertEqual(list(config.resources), ["gpu1", "gpu2", "gpu4", "gpu8x2"])
        self.assertEqual(config.resources["gpu4"].label, "4 GPUs")
        self.assertEqual(config.resources["gpu8x2"].label, "2 nodes / 16 GPUs")
        self.assertEqual(config.cases[-1].runner, "cusolvermp")

    def test_resource_names_are_not_hardcoded(self):
        text = (ROOT / "gpu-matrix.ini").read_text(encoding="utf-8")
        text = text.replace("resource.gpu1", "resource.single", 1)
        text = text.replace("resource = gpu1", "resource = single", 1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.ini"
            path.write_text(text, encoding="utf-8")
            self.assertIn("single", sai.load_config(path).resources)

    def test_unknown_resource_and_noncontiguous_case_fail(self):
        original = (ROOT / "gpu-matrix.ini").read_text(encoding="utf-8")
        for text in (
            original.replace("resource = gpu1", "resource = missing", 1),
            original.replace("[case.002]", "[case.999]", 1),
        ):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "matrix.ini"
                path.write_text(text, encoding="utf-8")
                with self.assertRaises(ValueError):
                    sai.load_config(path)


class CliTests(unittest.TestCase):
    def test_top_level_help_only_exposes_public_commands(self):
        help_text = subprocess.run(
            (sys.executable, str(ROOT / "sai.py"), "--help"),
            check=True, text=True, capture_output=True,
        ).stdout
        self.assertIn("{run,config}", help_text)
        self.assertIn("build and run the GPU matrix on SAI", help_text)
        self.assertNotIn("remote-run", help_text)
        self.assertNotIn("github-admit", help_text)

    def test_run_help_describes_every_option(self):
        help_text = subprocess.run(
            (sys.executable, str(ROOT / "sai.py"), "run", "--help"),
            check=True, text=True, capture_output=True,
        ).stdout
        for option in (
            "--ssh-config", "--target", "--project-root", "--source-repository",
            "--source-sha", "--namespace", "--run-id", "--run-attempt", "--artifacts",
        ):
            self.assertIn(option, help_text)
        self.assertIn("40-character commit SHA", help_text)
        self.assertIn("checkout's", help_text)

    def test_run_options_parse(self):
        arguments = [
            "run", "--ssh-config", "/tmp/ssh/config", "--target", "sai-ci",
            "--project-root", "/home/user/abacus_sai_gpu_ci",
            "--source-repository", "/tmp/abacus", "--source-sha", "a" * 40,
            "--namespace", "manual", "--run-id", "42", "--run-attempt", "1",
            "--artifacts", "/tmp/results",
        ]
        args = sai.parser().parse_args(arguments)
        self.assertEqual(args.command, "run")
        self.assertEqual(args.ssh_config, Path("/tmp/ssh/config"))
        self.assertEqual(args.source_repository, Path("/tmp/abacus"))
        self.assertEqual(args.artifacts, Path("/tmp/results"))
        self.assertEqual(args.source_sha, "a" * 40)

        for index in range(1, len(arguments), 2):
            with self.subTest(option=arguments[index]), \
                    mock.patch("sys.stderr", new_callable=io.StringIO), \
                    self.assertRaises(SystemExit):
                sai.parser().parse_args(arguments[:index] + arguments[index + 2:])


class TemplateTests(unittest.TestCase):
    def test_resource_template_is_complete_and_has_no_cpu_request(self):
        config = sai.load_config()
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "runs" / "manual" / "1-1"
            destination = run / "jobs" / "gpu4.sbatch"
            values = sai._job_values(
                config.resources["gpu4"], config, run, "gpu4",
                run / "results" / "gpu4-%A_%a.out",
            )
            values.update({
                "ARRAY": "0-3%4", "BUILD_JOB": "123",
                "MAPPING_ROOT": config.mapping_root,
                "DISABLE_NCCL_IB": "false", "MANIFEST": run / "jobs" / "gpu4.tsv",
            })
            sai._render(ROOT / "case.sbatch.in", destination, values)
            text = destination.read_text(encoding="utf-8")
            self.assertIn("#SBATCH --nodes=1", text)
            self.assertIn("#SBATCH --array=0-3%4", text)
            self.assertIn("#SBATCH --dependency=afterok:123", text)
            self.assertIn("SLURM_EXPORT_ENV=ALL", text)
            self.assertIn("OMPI_MCA_plm_slurm_args=--external-launcher", text)
            self.assertIn("PRTE_MCA_plm_slurm_args=--external-launcher", text)
            self.assertNotRegex(text, r"cpus-per-task|--mem(?:ory)?|--nodelist")
            self.assertNotRegex(text, r"@[A-Z_]+@")

    def test_modules_do_not_spell_dependency_paths(self):
        text = (ROOT / "modules.sh").read_text(encoding="utf-8")
        self.assertIn("module load abacus/", text)
        self.assertIn("LD_PRELOAD=${LD_PRELOAD:-}", text)
        self.assertIn("CMAKE_LIBRARY_PATH=${LIBRARY_PATH:-}", text)
        self.assertIn("CMAKE_INCLUDE_PATH=${CPATH:-}", text)
        self.assertNotRegex(text, r"CUSOLVERMP_PATH|CUBLASMP_PATH|NCCL_PATH|/lib/lib")


class SlurmTests(unittest.TestCase):
    def test_submit_and_accounting_require_each_array_task(self):
        responses = [
            mock.Mock(returncode=0, stdout="101\n", stderr=""),
            mock.Mock(returncode=0, stdout="", stderr=""),
            mock.Mock(
                returncode=0,
                stdout="101|COMPLETED|0:0\n101_0|COMPLETED|0:0\n101_1|FAILED|1:0\n",
                stderr="",
            ),
        ]
        with mock.patch("slurm.subprocess.run", side_effect=responses):
            client = slurm.Slurm(poll_seconds=0)
            job = client.submit(Path("case.sbatch"), array_count=2)
            states = client.wait((job,))
        self.assertEqual(states["101_0"], ("COMPLETED", "0:0"))
        self.assertEqual(states["101_1"], ("FAILED", "1:0"))

    def test_pass_requires_successful_slurm_accounting(self):
        config = sai.Config(
            "gpu", Path("/opt/sai_config/mps_mapping.d"), False, 1,
            sai.Resource("build", "flood-gpu", 1, 1, 1, 60),
            {"one": sai.Resource("one", "flood-gpu", 1, 1, 1, 60)},
            (sai.Case("suite", "case", "one", "autotest"),),
        )
        client = mock.Mock()
        client.submit.side_effect = ["100", "101"]
        client.wait.side_effect = [
            {"100": ("COMPLETED", "0:0")},
            {"101_0": ("FAILED", "1:0")},
        ]
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            (run / "jobs").mkdir()
            status = run / "results" / "status" / "suite__case.json"
            status.parent.mkdir(parents=True)
            status.write_text(json.dumps({"state": "PASS", "exit_code": 0}), encoding="utf-8")
            with mock.patch("sai.load_config", return_value=config), \
                    mock.patch("sai.Slurm", return_value=client), \
                    mock.patch("sai._render"):
                self.assertEqual(sai.remote_run(run), 1)
            row = json.loads((run / "results" / "result.json").read_text())["cases"][0]
            self.assertEqual(row["state"], "INFRA")
            self.assertEqual(row["slurm_exit_code"], "1:0")


class TransferTests(unittest.TestCase):
    @staticmethod
    def _git(root, *arguments):
        return subprocess.run(("git", *arguments), cwd=str(root), check=True, text=True, capture_output=True)

    def _bundle(self, repository, run, revision):
        bundle = run / "source.bundle.local"
        self._git(repository, "bundle", "create", str(bundle), revision)
        parts, checksum = sai._split_bundle(bundle, run / "input")
        bundle.unlink()
        self.assertEqual(len(parts), 8)
        return checksum

    def test_full_then_incremental_bundle(self):
        with tempfile.TemporaryDirectory() as directory:
            home = Path(directory)
            repository = home / "local"
            repository.mkdir()
            self._git(repository, "init")
            self._git(repository, "config", "user.email", "ci@example.invalid")
            self._git(repository, "config", "user.name", "CI")
            (repository / "value.txt").write_text("one\n", encoding="utf-8")
            self._git(repository, "add", ".")
            self._git(repository, "commit", "-m", "one")
            first = self._git(repository, "rev-parse", "HEAD").stdout.strip()

            project = home / "project"
            run1 = project / "runs" / "manual" / "1-1"
            (run1 / "control").mkdir(parents=True)
            with mock.patch("sai.Path.home", return_value=home):
                sai.remote_prepare(project, run1)
            checksum = self._bundle(repository, run1, "HEAD")
            sai.remote_receive(project, run1, first, checksum)
            self.assertEqual((run1 / "source" / "value.txt").read_text(), "one\n")

            (repository / "value.txt").write_text("two\n", encoding="utf-8")
            self._git(repository, "commit", "-am", "two")
            second = self._git(repository, "rev-parse", "HEAD").stdout.strip()
            run2 = project / "runs" / "manual" / "2-1"
            (run2 / "control").mkdir(parents=True)
            with mock.patch("sai.Path.home", return_value=home):
                sai.remote_prepare(project, run2)
            checksum = self._bundle(repository, run2, first + "..HEAD")
            sai.remote_receive(project, run2, second, checksum)
            self.assertEqual((run2 / "source" / "value.txt").read_text(), "two\n")

            run3 = project / "runs" / "manual" / "3-1"
            (run3 / "control").mkdir(parents=True)
            with mock.patch("sai.Path.home", return_value=home):
                sai.remote_prepare(project, run3)
            self.assertIsNone(sai._bundle_revision(repository, second, second))
            sai.remote_receive(project, run3, second, "-")
            self.assertEqual((run3 / "source" / "value.txt").read_text(), "two\n")

    def test_bundle_merge_rejects_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bundle = root / "bundle"
            bundle.write_bytes(bytes(range(256)) * 4)
            parts, checksum = sai._split_bundle(bundle, root)
            parts[0].write_bytes(b"corrupt")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                sai._assemble_bundle(root, checksum)

    def test_prepare_rejects_reused_run(self):
        with tempfile.TemporaryDirectory() as directory:
            home = Path(directory)
            project = home / "project"
            run = project / "runs" / "manual" / "1-1"
            (run / "control").mkdir(parents=True)
            with mock.patch("sai.Path.home", return_value=home):
                sai.remote_prepare(project, run)
                with self.assertRaisesRegex(ValueError, "already exists"):
                    sai.remote_prepare(project, run)

    def test_transfer_retry_is_bounded(self):
        completed = mock.Mock(returncode=0, stdout="done", stderr="")
        failed = mock.Mock(returncode=1, stdout="", stderr="disconnected")
        with mock.patch("sai.subprocess.run", side_effect=[failed, failed, completed]) as run, \
                mock.patch("sai.time.sleep"):
            self.assertEqual(sai._retry(("rsync", "source", "target")).stdout, "done")
        self.assertEqual(run.call_count, 3)

    def test_download_retry_replaces_partial_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.tar.gz"

            def download(_command, _cwd=None, stdout=None):
                stdout.write(b"partial" if download.calls == 0 else b"complete")
                download.calls += 1
                if download.calls == 1:
                    raise RuntimeError("connection closed")

            download.calls = 0
            with mock.patch("sai._command", side_effect=download), mock.patch("sai.time.sleep"):
                sai._retry_download(("ssh", "collect"), path)
            self.assertEqual(path.read_bytes(), b"complete")

    def test_run_records_cleanup_metadata_before_remote_work(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_sha = "a" * 40
            args = argparse.Namespace(
                source_repository=root,
                project_root="/project",
                target="sai-ci",
                namespace="manual",
                run_id="1",
                run_attempt="1",
                source_sha=source_sha,
                artifacts=root / "artifacts",
                ssh_config=root / "ssh-config",
            )
            revision = mock.Mock(stdout=source_sha + "\n")
            with mock.patch("sai._command", return_value=revision), \
                    mock.patch("sai._retry", side_effect=RuntimeError("disconnected")), \
                    self.assertRaisesRegex(RuntimeError, "disconnected"):
                sai.run(args)
            metadata = json.loads((args.artifacts / "run.json").read_text())
            self.assertEqual(metadata, {
                "project_root": "/project",
                "run_root": "/project/runs/manual/1-1",
                "source_sha": source_sha,
            })

    def test_result_archive_rejects_traversal_and_links(self):
        for name, link in (("../sai-ssh/id_ed25519", None), ("results/key", "../key")):
            payload = io.BytesIO()
            with tarfile.open(fileobj=payload, mode="w:gz") as archive:
                member = tarfile.TarInfo(name)
                if link is None:
                    member.size = 3
                    archive.addfile(member, io.BytesIO(b"key"))
                else:
                    member.type = tarfile.SYMTYPE
                    member.linkname = link
                    archive.addfile(member)
            payload.seek(0)
            with tempfile.TemporaryDirectory() as directory, self.assertRaises(ValueError):
                sai._extract_results(payload, Path(directory))

    def test_result_archive_extracts_regular_results(self):
        payload = io.BytesIO()
        with tarfile.open(fileobj=payload, mode="w:gz") as archive:
            member = tarfile.TarInfo("results/result.json")
            member.size = 3
            archive.addfile(member, io.BytesIO(b"{}\n"))
        payload.seek(0)
        with tempfile.TemporaryDirectory() as directory:
            sai._extract_results(payload, Path(directory))
            self.assertEqual((Path(directory) / "results" / "result.json").read_text(), "{}\n")

    def test_result_archive_rejects_oversized_content(self):
        payload = io.BytesIO()
        with tarfile.open(fileobj=payload, mode="w:gz") as archive:
            member = tarfile.TarInfo("results/large")
            member.size = 3
            archive.addfile(member, io.BytesIO(b"abc"))
        payload.seek(0)
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch("sai.MAX_RESULT_BYTES", 2), \
                self.assertRaisesRegex(ValueError, "too large"):
            sai._extract_results(payload, Path(directory))

    def test_result_archive_limits_directory_members(self):
        payload = io.BytesIO()
        with tarfile.open(fileobj=payload, mode="w:gz") as archive:
            for name in ("results/one", "results/two"):
                member = tarfile.TarInfo(name)
                member.type = tarfile.DIRTYPE
                archive.addfile(member)
        payload.seek(0)
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch("sai.MAX_RESULT_MEMBERS", 1), \
                self.assertRaisesRegex(ValueError, "too large"):
            sai._extract_results(payload, Path(directory))


class ResultTests(unittest.TestCase):
    def test_report_publishes_dynamic_components(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = valid_result()
            path = root / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            args = argparse.Namespace(result=path, output=root / "output", summary=root / "summary")
            sai.report(args)
            output = args.output.read_text(encoding="utf-8")
            self.assertIn("available=true", output)
            self.assertIn('"name":"gpu8x2"', output)
            self.assertTrue(args.summary.read_text().startswith("# SAI GPU result\n"))

    def test_report_rejects_untrusted_counts(self):
        invalid = (
            {"passed": "x[$(printf ARITH_EXEC >&2)0]", "failed": 0, "infrastructure": 0, "total": 1},
            {"passed": True, "failed": 0, "infrastructure": 0, "total": 1},
            {"passed": -1, "failed": 1, "infrastructure": 0, "total": 0},
            {"passed": 1, "failed": 1, "infrastructure": 0, "total": 1},
        )
        for counts in invalid:
            with self.subTest(counts=counts), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                path = root / "result.json"
                result = valid_result()
                result.update(counts)
                path.write_text(json.dumps(result), encoding="utf-8")
                args = argparse.Namespace(result=path, output=root / "output", summary=None)
                with self.assertRaisesRegex(ValueError, "result counts"):
                    sai.report(args)
                self.assertFalse(args.output.exists())

    def test_report_rejects_wrong_matrix_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = valid_result()
            result["cases"][0]["case_id"] = "invented/case"
            path = root / "result.json"
            path.write_text(json.dumps(result), encoding="utf-8")
            args = argparse.Namespace(result=path, output=root / "output", summary=None)
            with self.assertRaisesRegex(ValueError, "result case"):
                sai.report(args)

    def test_mpi_startup_failure_requires_complete_signature(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log"
            path.write_bytes(b"PMIX_ERR_FILE_OPEN_FAILURE MPI_Init_thread PMIx_Init failed")
            self.assertTrue(sai._mpi_startup_failure(path))
            path.write_bytes(b"PMIX_ERR_FILE_OPEN_FAILURE")
            self.assertFalse(sai._mpi_startup_failure(path))
            path.write_bytes(b"srun returned non-zero exit status (512) from launching the per-node daemon")
            self.assertTrue(sai._mpi_startup_failure(path))
            path.write_bytes(b"srun returned non-zero exit status (512)")
            self.assertFalse(sai._mpi_startup_failure(path))


class PolicyTests(unittest.TestCase):
    def test_workflow_uses_trusted_control_and_protected_environment(self):
        text = (ROOT.parents[1] / ".github" / "workflows" / "sai-gpu.yml").read_text(encoding="utf-8")
        self.assertIn("ref: ${{ env.CONTROL_SHA }}", text)
        self.assertIn("sai-ssh-scheduled", text)
        self.assertIn("sai-ssh-manual", text)
        self.assertIn("/abacus-ci sai-gpu", text)
        self.assertIn("sai.py run", text)
        self.assertIn("pull-requests: write", text)
        self.assertNotIn("cpus-per-task", text)


if __name__ == "__main__":
    unittest.main()
