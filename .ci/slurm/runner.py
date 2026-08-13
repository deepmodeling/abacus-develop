#!/usr/bin/env python3
"""Run the ABACUS GPU matrix on a remote Slurm cluster."""

import argparse
import configparser
import fcntl
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from pathlib import PurePosixPath
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import urlsplit

from slurm import Slurm, SlurmError


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parents[1]
NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
SHA = re.compile(r"[0-9a-f]{40}\Z")
PMIX = re.compile(br"PMIX_ERR_(?:FILE_OPEN_FAILURE|OUT_OF_RESOURCE)")
SRUN_DAEMON = re.compile(
    br"srun returned non-zero exit status \([1-9][0-9]*\) from launching\s+the per-node daemon"
)
MAX_RESULT_MEMBERS = 10000
MAX_RESULT_BYTES = 2 * 1024**3
MAX_REPORT_BYTES = 1024**2
BUNDLE_PARTS = 8
RECENT_CACHE_LIMIT = 3
DAILY_CACHE_LIMIT = 2
CACHE_RESERVATION_SECONDS = 6 * 3600


@dataclass(frozen=True)
class Resource:
    name: str
    qos: str
    nodes: int
    tasks_per_node: int
    gpus_per_node: int
    time_seconds: int
    parallelism: int = 1

    @property
    def tasks(self) -> int:
        return self.nodes * self.tasks_per_node

    @property
    def label(self) -> str:
        total = self.nodes * self.gpus_per_node
        if self.nodes == 1:
            return "{} GPU{}".format(total, "" if total == 1 else "s")
        return "{} nodes / {} GPUs".format(self.nodes, total)


@dataclass(frozen=True)
class Case:
    suite: str
    name: str
    resource: str
    runner: str

    @property
    def case_id(self) -> str:
        return self.suite + "/" + self.name


@dataclass(frozen=True)
class Remote:
    host: str
    port: int
    user: str
    project_root: str
    allowed_project_roots: Tuple[PurePosixPath, ...]


@dataclass(frozen=True)
class Site:
    name: str
    url: str
    acknowledgement: str


@dataclass(frozen=True)
class Config:
    site: Site
    remote: Remote
    partition: str
    mapping_root: Path
    disable_nccl_ib: bool
    poll_seconds: int
    build: Resource
    resources: Mapping[str, Resource]
    cases: Tuple[Case, ...]


def _integer(section: Mapping[str, str], key: str, low: int, high: int) -> int:
    try:
        value = int(section[key])
    except (KeyError, ValueError) as error:
        raise ValueError("invalid {}".format(key)) from error
    if not low <= value <= high:
        raise ValueError("{} is outside {}..{}".format(key, low, high))
    return value


def _resource(name: str, section: Mapping[str, str], array: bool) -> Resource:
    expected = {
        "qos", "nodes", "tasks_per_node", "gpus_per_node", "time_seconds",
    } | ({"parallelism"} if array else set())
    if set(section) != expected or not NAME.fullmatch(name):
        raise ValueError("invalid resource {}".format(name))
    qos = section["qos"]
    if not NAME.fullmatch(qos):
        raise ValueError("invalid QoS")
    profile = Resource(
        name, qos, _integer(section, "nodes", 1, 2),
        _integer(section, "tasks_per_node", 1, 8),
        _integer(section, "gpus_per_node", 1, 8),
        _integer(section, "time_seconds", 1, 3600),
        _integer(section, "parallelism", 1, 32) if array else 1,
    )
    if profile.tasks_per_node != profile.gpus_per_node or profile.tasks > 16:
        raise ValueError("{} must use one rank per GPU and at most 16 GPUs".format(name))
    return profile


def load_config(path: Path = ROOT / "config.ini") -> Config:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    with Path(path).open(encoding="utf-8") as stream:
        parser.read_file(stream)
    resources = [name for name in parser.sections() if name.startswith("resource.")]
    cases = [name for name in parser.sections() if name.startswith("case.")]
    expected_cases = ["case.{:03d}".format(index) for index in range(1, len(cases) + 1)]
    if not resources or cases != expected_cases:
        raise ValueError("resources and contiguous case sections are required")
    known = {"site", "remote", "cluster", "build", *resources, *cases}
    if set(parser.sections()) != known or set(parser["site"]) != {
        "name", "url", "acknowledgement",
    } or set(parser["remote"]) != {
        "host", "port", "user", "project_root", "allowed_project_roots",
    } or set(parser["cluster"]) != {
        "partition", "mapping_root", "disable_nccl_ib", "poll_seconds",
    }:
        raise ValueError("unexpected configuration section or key")
    site = parser["site"]
    site_url = urlsplit(site["url"])
    if not site["name"].strip() or not site["acknowledgement"].strip() or \
            "\n" in site["name"] or "\n" in site["acknowledgement"] or \
            site_url.scheme != "https" or not site_url.netloc:
        raise ValueError("invalid site configuration")
    remote = parser["remote"]
    project_root = remote["project_root"]
    project = PurePosixPath(project_root)
    if project_root == "~":
        project_parts: Tuple[str, ...] = ()
    elif project_root.startswith("~/"):
        project_parts = project.parts[1:]
    elif project.is_absolute():
        project_parts = project.parts[1:]
    else:
        raise ValueError("invalid remote project root")
    allowed_values = [value.strip() for value in remote["allowed_project_roots"].split(",")]
    allowed_roots = tuple(PurePosixPath(value) for value in allowed_values)
    if not NAME.fullmatch(remote["host"]) or not NAME.fullmatch(remote["user"]) or \
            any(not NAME.fullmatch(part) for part in project_parts) or not allowed_roots or \
            len(set(allowed_roots)) != len(allowed_roots) or any(
                not root.is_absolute() or len(root.parts) < 2 or
                any(not NAME.fullmatch(part) for part in root.parts[1:])
                for root in allowed_roots
            ):
        raise ValueError("invalid remote configuration")
    partition = parser["cluster"]["partition"]
    mapping = Path(parser["cluster"]["mapping_root"])
    disable = parser["cluster"]["disable_nccl_ib"].lower()
    if not NAME.fullmatch(partition) or not mapping.is_absolute() or disable not in ("true", "false"):
        raise ValueError("invalid cluster configuration")
    profiles: "OrderedDict[str, Resource]" = OrderedDict()
    for section_name in resources:
        name = section_name[len("resource."):]
        profiles[name] = _resource(name, parser[section_name], True)
    matrix: List[Case] = []
    for section_name in cases:
        section = parser[section_name]
        if set(section) != {"suite", "name", "resource", "runner"}:
            raise ValueError("invalid {}".format(section_name))
        case = Case(section["suite"], section["name"], section["resource"], section["runner"])
        if not all(NAME.fullmatch(value) for value in (case.suite, case.name, case.resource)):
            raise ValueError("invalid case name")
        if case.resource not in profiles or case.runner not in (
            "autotest", "autotest_gpu", "cusolvermp",
        ):
            raise ValueError("invalid case resource or runner")
        matrix.append(case)
    if len({case.case_id for case in matrix}) != len(matrix):
        raise ValueError("duplicate case")
    if any(not any(case.resource == name for case in matrix) for name in profiles):
        raise ValueError("every resource needs a case")
    return Config(
        Site(site["name"], site["url"], site["acknowledgement"]),
        Remote(
            remote["host"], _integer(remote, "port", 1, 65535),
            remote["user"], project_root, allowed_roots,
        ),
        partition, mapping, disable == "true",
        _integer(parser["cluster"], "poll_seconds", 1, 300),
        _resource("build", parser["build"], False), profiles, tuple(matrix),
    )


def _atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".{}.{}.tmp".format(path.name, os.getpid()))
    text = value if isinstance(value, str) else json.dumps(value, indent=2, sort_keys=True) + "\n"
    temporary.write_text(text, encoding="utf-8")
    os.replace(str(temporary), str(path))


def _below_project_root(project: Path, roots: Sequence[Path]) -> bool:
    return bool(roots) and all(root.is_absolute() and len(root.parts) > 1 for root in roots) and \
        any(root in project.parents for root in roots)


def _time(seconds: int) -> str:
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02d}:{:02d}:{:02d}".format(hours, minutes, seconds)


def _render(template: Path, destination: Path, values: Mapping[str, Any]) -> None:
    text = template.read_text(encoding="utf-8")
    for name, value in values.items():
        text = text.replace("@{}@".format(name), str(value))
    if re.search(r"@[A-Z_]+@", text):
        raise ValueError("unresolved Slurm template value")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    destination.chmod(0o755)


def _job_values(
    profile: Resource, config: Config, run: Path, label: str, output: Path,
) -> Dict[str, Any]:
    return {
        "JOB_NAME": "abacus-{}-{}".format(run.name, label),
        "PARTITION": config.partition, "QOS": profile.qos,
        "NODES": profile.nodes, "TASKS": profile.tasks,
        "TASKS_PER_NODE": profile.tasks_per_node,
        "GPUS_PER_NODE": profile.gpus_per_node,
        "TIME": _time(profile.time_seconds), "OUTPUT": output,
        "HOME": Path.home(), "CONTROL": run / "control",
        "SOURCE": run / "source", "BUILD": run / "build",
        "INSTALL": run / "install", "RESULTS": run / "results",
    }


def _component(name: str, label: str, state: str, job: str = "", slurm: str = "", code: str = "") -> Dict[str, str]:
    return {
        "name": name, "label": label, "state": state, "job_id": job,
        "slurm_state": slurm, "exit_code": code,
    }


def _folder_components(result: Mapping[str, Any]) -> List[Dict[str, str]]:
    components = [dict(result["components"][0])]
    suites: Dict[str, List[Mapping[str, Any]]] = {}
    for row in result["cases"]:
        suite = row["case_id"].split("/", 1)[0]
        suites.setdefault(suite, []).append(row)
    for suite in sorted(suites):
        rows = suites[suite]
        states = [row["state"] for row in rows]
        state = "PASS" if all(item == "PASS" for item in states) else (
            "FAIL" if any(item in ("FAIL", "TIMEOUT") for item in states) else "INFRA"
        )
        jobs = list(dict.fromkeys(
            row["job_id"].split("_", 1)[0] for row in rows if row["job_id"]
        ))
        components.append(_component(
            suite, "tests/" + suite, state, ", ".join(jobs),
        ))
    return components


def _result_row(case: Case, state: str, **values: Any) -> Dict[str, Any]:
    row = {
        "case_id": case.case_id, "resource": case.resource,
        "runner": case.runner, "state": state, "exit_code": None,
        "slurm_state": "", "slurm_exit_code": "", "job_id": "",
        "elapsed_seconds": 0,
    }
    row.update(values)
    return row


def _save_result(run: Path, components: Sequence[Mapping[str, Any]], rows: Sequence[Mapping[str, Any]]) -> int:
    passed = sum(row["state"] == "PASS" for row in rows)
    failed = sum(row["state"] in ("FAIL", "TIMEOUT") for row in rows)
    result = {
        "protocol": 1, "total": len(rows), "passed": passed, "failed": failed,
        "infrastructure": len(rows) - passed - failed,
        "components": list(components), "cases": list(rows),
    }
    root = run / "results"
    _atomic(root / "result.json", result)
    _atomic(root / "summary.md", _result_markdown(result))
    return 0 if rows and passed == len(rows) else 1


def _site_credit() -> str:
    site = load_config().site
    return "{} [{}]({}).".format(site.acknowledgement, site.name, site.url)


def _result_markdown(result: Mapping[str, Any]) -> str:
    components = _folder_components(result)
    lines = [
        "# GPU validation result", "",
        "Passed: **{}**; failed: **{}**; infrastructure: **{}**".format(
            result["passed"], result["failed"], result["infrastructure"]
        ), "", "| Component | State | Slurm jobs |", "| --- | --- | --- |",
    ]
    lines.extend("| {} | {} | {} |".format(item["label"], item["state"], item["job_id"]) for item in components)
    lines.extend(("", "| Case | Resource | State | Duration | Slurm job |", "| --- | --- | --- | --- | --- |"))
    lines.extend("| {} | {} | {} | {} | {} |".format(
        row["case_id"], row["resource"], row["state"],
        _time(row["elapsed_seconds"]), row["job_id"],
    ) for row in sorted(result["cases"], key=lambda item: item["case_id"]))
    lines.extend(("", _site_credit()))
    return "\n".join(lines) + "\n"


def remote_run(run: Path) -> int:
    run = run.resolve()
    config = load_config(run / "control" / "config.ini")
    results = run / "results"
    scripts = run / "jobs"
    results.mkdir(parents=True, exist_ok=True)
    slurm = Slurm(config.poll_seconds)
    components: List[Mapping[str, Any]] = []
    rows: List[Mapping[str, Any]] = []
    returncode = 2
    try:
        build_script = scripts / "build.sbatch"
        _render(
            run / "control" / "build.sbatch.in", build_script,
            _job_values(config.build, config, run, "build", results / "build-%j.out"),
        )
        build_job = slurm.submit(build_script)
        build_jobs = (("Compile", build_job),)
        build_progress: Dict[str, Tuple[int, int, int]] = {}
        build_state, build_exit = slurm.wait(
            (build_job,), lambda states: _print_progress(build_jobs, states, build_progress),
        )[build_job]
        build_ok = (build_state, build_exit) == ("COMPLETED", "0:0")
        components.append(_component(
            "build", "Compile", "PASS" if build_ok else "FAIL",
            build_job, build_state, build_exit,
        ))
        if not build_ok:
            components.extend(_component(name, profile.label, "SKIPPED") for name, profile in config.resources.items())
            rows = [
                _result_row(
                    case, "INFRA", job_id=build_job, slurm_state=build_state,
                    slurm_exit_code=build_exit,
                )
                for resource in config.resources
                for case in config.cases if case.resource == resource
            ]
            returncode = _save_result(run, components, rows)
            return returncode

        jobs: Dict[str, str] = {}
        grouped: Dict[str, List[Case]] = {}
        for name, profile in config.resources.items():
            grouped[name] = [case for case in config.cases if case.resource == name]
            manifest = scripts / (name + ".tsv")
            manifest.write_text("\n".join("\t".join((case.case_id, case.suite, case.name, case.resource, case.runner)) for case in grouped[name]) + "\n", encoding="utf-8")
            values = _job_values(profile, config, run, name, results / (name + "-%A_%a.out"))
            values.update({
                "ARRAY": "0-{}%{}".format(len(grouped[name]) - 1, min(profile.parallelism, len(grouped[name]))),
                "BUILD_JOB": build_job, "MAPPING_ROOT": config.mapping_root,
                "DISABLE_NCCL_IB": str(config.disable_nccl_ib).lower(),
                "MANIFEST": manifest,
            })
            script = scripts / (name + ".sbatch")
            _render(run / "control" / "case.sbatch.in", script, values)
            jobs[name] = slurm.submit(script, len(grouped[name]))
        test_jobs = tuple((config.resources[name].label, job) for name, job in jobs.items())
        test_progress: Dict[str, Tuple[int, int, int]] = {}
        accounting = slurm.wait(
            tuple(jobs.values()), lambda states: _print_progress(test_jobs, states, test_progress),
        )

        for name, profile in config.resources.items():
            group_states = []
            for index, case in enumerate(grouped[name]):
                job = "{}_{}".format(jobs[name], index)
                slurm_state, slurm_exit = accounting[job]
                status = results / "status" / (case.case_id.replace("/", "__") + ".json")
                try:
                    data = json.loads(status.read_text(encoding="utf-8"))
                    state = data["state"]
                    if state not in ("PASS", "FAIL", "TIMEOUT", "INFRA"):
                        raise ValueError
                except (OSError, ValueError, KeyError, json.JSONDecodeError):
                    data, state = {}, "INFRA"
                if state == "PASS" and (slurm_state, slurm_exit) != ("COMPLETED", "0:0"):
                    state = "INFRA"
                rows.append(_result_row(
                    case, state, exit_code=data.get("exit_code"),
                    elapsed_seconds=data.get("elapsed_seconds", 0),
                    slurm_state=slurm_state, slurm_exit_code=slurm_exit,
                    job_id=job,
                ))
                group_states.append(state)
            state = "PASS" if all(item == "PASS" for item in group_states) else (
                "FAIL" if any(item in ("FAIL", "TIMEOUT") for item in group_states) else "INFRA"
            )
            components.append(_component(name, profile.label, state, jobs[name]))
        returncode = _save_result(run, components, rows)
        return returncode
    except Exception as error:
        slurm.cancel()
        _atomic(results / "coordinator-error.txt", str(error) + "\n")
        print("gpu-ci: {}".format(error), file=sys.stderr, flush=True)
        return returncode
    finally:
        _atomic(results / "done.json", {"returncode": returncode})


def _stream(command: Sequence[str], cwd: Path, log: Path) -> int:
    with log.open("wb") as output:
        process = subprocess.Popen(command, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.stdout is not None
        for block in iter(lambda: process.stdout.read(65536), b""):
            output.write(block)
            output.flush()
            sys.stdout.buffer.write(block)
            sys.stdout.buffer.flush()
        return process.wait()


def _mpi_startup_failure(log: Path) -> bool:
    data = log.read_bytes()
    return (
        bool(PMIX.search(data)) and b"MPI_Init_thread" in data
    ) or bool(SRUN_DAEMON.search(data))


def _force_gpu_inputs(case: Path, artifacts: Optional[Path] = None) -> None:
    inputs = sorted(path for path in case.rglob("INPUT") if path.is_file())
    if not inputs:
        raise ValueError("GPU autotest case has no INPUT")
    for path in inputs:
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines(keepends=True)
        active = [
            index for index, line in enumerate(lines)
            if re.match(r"^\s*device(?:\s*=|\s+)", line) and not line.lstrip().startswith("#")
        ]
        if len(active) > 1:
            raise ValueError("GPU autotest INPUT has duplicate device entries: {}".format(path))
        if active:
            index = active[0]
            ending = "\r\n" if lines[index].endswith("\r\n") else "\n" if lines[index].endswith("\n") else ""
            body = lines[index][:-len(ending)] if ending else lines[index]
            match = re.match(
                r"^(\s*)device(?:\s*=|\s+)\s*\S+(\s*(?:#.*)?)$", body,
            )
            if not match:
                raise ValueError("invalid device entry in GPU autotest INPUT: {}".format(path))
            lines[index] = "{}device gpu{}{}".format(match.group(1), match.group(2), ending)
        else:
            if text and not text.endswith(("\n", "\r")):
                lines.append("\n")
            lines.append("device gpu\n")
        path.write_text("".join(lines), encoding="utf-8")
        if artifacts is not None:
            destination = artifacts / "effective-inputs" / path.relative_to(case)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(str(path), str(destination))


def worker(source: Path, control: Path, install: Path, results: Path, manifest: Path) -> int:
    task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    fields = manifest.read_text(encoding="utf-8").splitlines()[task_id].split("\t")
    if len(fields) != 5:
        raise ValueError("invalid case manifest")
    case_id, suite, name, resource, runner = fields
    artifacts = results / "cases" / case_id.replace("/", "__")
    artifacts.mkdir(parents=True, exist_ok=True)
    source = source.resolve()
    work = results.parent / "work" / case_id.replace("/", "__")
    tests = work / "tests"
    case = tests / suite / name
    shutil.rmtree(str(work), ignore_errors=True)
    (tests / suite).mkdir(parents=True)
    os.symlink(str(source / "tests" / "integrate"), str(tests / "integrate"))
    os.symlink(str(source / "tests" / "PP_ORB"), str(tests / "PP_ORB"))
    shutil.copytree(str(source / "tests" / suite / name), str(case))
    if runner == "autotest_gpu":
        _force_gpu_inputs(case, artifacts)
    launcher = work / "launcher"
    launcher.mkdir()
    os.symlink(str(control / "mpirun_with_mapping.sh"), str(launcher / "mpirun"))
    os.environ["PATH"] = str(launcher) + os.pathsep + os.environ["PATH"]
    abacus = install / "bin" / "abacus"
    listing = subprocess.run(("ldd", str(abacus)), text=True, capture_output=True)
    (artifacts / "ldd.txt").write_text(listing.stdout + listing.stderr, encoding="utf-8")
    if listing.returncode or "not found" in listing.stdout:
        raise RuntimeError("ABACUS has unresolved runtime libraries")
    start = time.time()
    returncode = 2
    final_startup_failure = False
    try:
        if runner in ("autotest", "autotest_gpu"):
            cases_file = case.parent / "CASES.task.txt"
            cases_file.write_text(name + "\n", encoding="utf-8")
            command = (
                "timeout", "--signal=TERM", "--kill-after=30s", "10m",
                "bash", str(tests / "integrate" / "Autotest.sh"),
                "-a", str(abacus), "-n", os.environ["SLURM_NTASKS"],
                "-o", "1", "-f", cases_file.name, "-r", "^{}$".format(name),
            )
            for attempt in (1, 2):
                shutil.rmtree(str(case / "OUT.autotest"), ignore_errors=True)
                for filename in ("log.txt", "result.out", "warning.log"):
                    try:
                        (case / filename).unlink()
                    except FileNotFoundError:
                        pass
                log = artifacts / "attempt-{}.log".format(attempt)
                returncode = _stream(command, case.parent, log)
                final_startup_failure = returncode != 0 and _mpi_startup_failure(log)
                if returncode == 0 or attempt == 2 or not final_startup_failure:
                    break
                time.sleep(10)
        else:
            command = (
                "timeout", "--signal=TERM", "--kill-after=30s", "35m",
                "mpirun", "-np", os.environ["SLURM_NTASKS"], str(abacus),
            )
            for attempt in (1, 2):
                log = artifacts / ("cusolvermp.log" if attempt == 1 else "cusolvermp-retry.log")
                returncode = _stream(command, case, log)
                final_startup_failure = returncode != 0 and _mpi_startup_failure(log)
                if returncode == 0 or attempt == 2 or not final_startup_failure:
                    break
                time.sleep(10)
        for filename in ("log.txt", "result.out", "warning.log"):
            if (case / filename).is_file():
                shutil.copy2(str(case / filename), str(artifacts / filename))
        state = "PASS" if returncode == 0 else (
            "TIMEOUT" if returncode in (124, 137, 143) else
            "INFRA" if final_startup_failure else "FAIL"
        )
        _atomic(results / "status" / (case_id.replace("/", "__") + ".json"), {
            "case_id": case_id, "resource": resource, "runner": runner,
            "state": state, "exit_code": returncode,
            "elapsed_seconds": int(time.time() - start),
        })
        return returncode
    except Exception:
        _atomic(results / "status" / (case_id.replace("/", "__") + ".json"), {
            "case_id": case_id, "resource": resource, "runner": runner,
            "state": "INFRA", "exit_code": 2,
            "elapsed_seconds": int(time.time() - start),
        })
        raise


def _command(command: Sequence[str], cwd: Optional[Path] = None, stdout: Any = subprocess.PIPE) -> subprocess.CompletedProcess:
    result = subprocess.run(command, cwd=str(cwd) if cwd else None, text=True, stdout=stdout, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "command failed: {}".format(command[0]))
    return result


def _retry(command: Sequence[str], cwd: Optional[Path] = None, stdout: Any = subprocess.PIPE) -> subprocess.CompletedProcess:
    error: Optional[Exception] = None
    for attempt in range(3):
        try:
            return _command(command, cwd, stdout)
        except RuntimeError as caught:
            error = caught
            if attempt < 2:
                time.sleep(5 * (attempt + 1))
    raise RuntimeError(str(error))


def _retry_download(command: Sequence[str], destination: Path) -> None:
    error: Optional[Exception] = None
    for attempt in range(3):
        try:
            with destination.open("wb") as output:
                _command(command, stdout=output)
            return
        except RuntimeError as caught:
            error = caught
            if attempt < 2:
                time.sleep(5 * (attempt + 1))
    raise RuntimeError(str(error))


def _cache(project: Path) -> Path:
    return project / "cache" / "repository.git"


def _run_relative(project: Path, run: Path) -> Path:
    try:
        relative = run.resolve().relative_to((project.resolve() / "runs").resolve())
    except ValueError as error:
        raise ValueError("invalid run directory") from error
    if len(relative.parts) != 2 or any(not NAME.fullmatch(part) for part in relative.parts):
        raise ValueError("invalid run directory")
    return relative


def _cache_refs(cache: Path, prefix: str = "refs/ci") -> List[Tuple[str, str]]:
    lines = _command((
        "git", "--git-dir", str(cache), "for-each-ref",
        "--format=%(refname) %(objectname)", prefix,
    )).stdout.splitlines()
    refs = []
    for line in lines:
        try:
            name, value = line.split()
        except ValueError as error:
            raise ValueError("invalid source cache ref") from error
        if not SHA.fullmatch(value):
            raise ValueError("invalid cached source SHA")
        refs.append((name, value))
    return refs


def _retained_refs(cache: Path) -> List[Tuple[str, str]]:
    refs = []
    for ref, value in _cache_refs(cache):
        if ref.startswith("refs/ci/active/"):
            continue
        if not re.fullmatch(
            r"refs/ci/(?:daily/(?:day/\d{4}-\d{2}-\d{2}|\d+)|"
            r"weekly/\d{4}-\d{2}/\d{4}-W\d{2}|monthly/\d{4}-\d{2}|"
            r"recent/\d+|[0-9a-f]{40})",
            ref,
        ):
            raise ValueError("invalid source cache ref")
        refs.append((ref, value))
    return refs


def _cached_commits(cache: Path, retained: Sequence[Tuple[str, str]]) -> List[str]:
    refs = [name for name, _ in retained]
    if not refs:
        return []
    commits = _command(("git", "--git-dir", str(cache), "rev-list", *refs)).stdout.splitlines()
    if any(not SHA.fullmatch(value) for value in commits):
        raise ValueError("invalid cached source SHA")
    return sorted(set(commits))


def _update_cache_refs(cache: Path, updates: Mapping[str, str], deletes: Sequence[str]) -> None:
    if not updates and not deletes:
        return
    commands = ["start"]
    commands.extend("update {} {}".format(ref, value) for ref, value in updates.items())
    commands.extend("delete {}".format(ref) for ref in deletes)
    commands.extend(("prepare", "commit"))
    result = subprocess.run(
        ("git", "--git-dir", str(cache), "update-ref", "--stdin"),
        input="\n".join(commands) + "\n", text=True, capture_output=True,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "unable to update source cache refs")


def _reservation_token(run: Path) -> str:
    return hashlib.sha256(str(run.resolve()).encode("utf-8")).hexdigest()


def _active_refs(cache: Path) -> List[Tuple[str, int, str]]:
    active = []
    for ref, _ in _cache_refs(cache, "refs/ci/active"):
        match = re.fullmatch(r"refs/ci/active/(\d+)/([0-9a-f]{64})/\d+", ref)
        if not match:
            raise ValueError("invalid source cache reservation")
        active.append((ref, int(match.group(1)), match.group(2)))
    return active


def _cleanup_reservations(cache: Path) -> None:
    now = time.time()
    _update_cache_refs(cache, {}, [ref for ref, expires, _ in _active_refs(cache) if expires <= now])


def _reserve_cache(cache: Path, run: Path, retained: Sequence[Tuple[str, str]]) -> None:
    if not retained:
        return
    prefix = "refs/ci/active/{}/{}".format(
        int(time.time() + CACHE_RESERVATION_SECONDS), _reservation_token(run),
    )
    _update_cache_refs(cache, {
        "{}/{}".format(prefix, index): value
        for index, (_, value) in enumerate(retained)
    }, ())


def _release_cache(cache: Path, run: Path) -> None:
    token = _reservation_token(run)
    refs = [ref for ref, _, owner in _active_refs(cache) if owner == token]
    _update_cache_refs(cache, {}, refs)


def _rotate_cache_refs(cache: Path, bucket: str, source_sha: str) -> None:
    days: Dict[str, str] = {}
    weeks: Dict[Tuple[str, str], str] = {}
    months: Dict[str, str] = {}
    recent: List[Tuple[int, str]] = []
    legacy = []
    retained = _retained_refs(cache)
    for ref, value in retained:
        daily_match = re.fullmatch(r"refs/ci/daily/day/(\d{4}-\d{2}-\d{2})", ref)
        weekly_match = re.fullmatch(
            r"refs/ci/weekly/(\d{4}-\d{2})/(\d{4}-W\d{2})", ref,
        )
        monthly_match = re.fullmatch(r"refs/ci/monthly/(\d{4}-\d{2})", ref)
        recent_match = re.fullmatch(r"refs/ci/recent/(\d+)", ref)
        if daily_match:
            days[daily_match.group(1)] = value
        elif weekly_match:
            weeks[(weekly_match.group(1), weekly_match.group(2))] = value
        elif monthly_match:
            months[monthly_match.group(1)] = value
        elif recent_match:
            recent.append((int(recent_match.group(1)), value))
        elif ref == "refs/ci/" + source_sha:
            continue
        elif re.fullmatch(r"refs/ci/(?:daily/\d+|[0-9a-f]{40})", ref):
            legacy.append(value)
        else:
            raise ValueError("invalid source cache ref")

    stamp = time.gmtime()
    today = date(stamp.tm_year, stamp.tm_mon, stamp.tm_mday)
    month = today.strftime("%Y-%m")
    if bucket == "daily":
        iso = today.isocalendar()
        week = "{}-W{:02d}".format(iso.year, iso.week)
        days[today.isoformat()] = source_sha
        weeks.setdefault((month, week), source_sha)
        months.setdefault(month, source_sha)
    else:
        recent.append((-1, source_sha))

    desired: Dict[str, str] = {}
    for name, value in sorted(months.items()):
        desired["refs/ci/monthly/{}".format(name)] = value
    for (week_month, week), value in sorted(weeks.items()):
        if week_month == month:
            desired["refs/ci/weekly/{}/{}".format(month, week)] = value
    for day in sorted(days, reverse=True)[:DAILY_CACHE_LIMIT]:
        desired["refs/ci/daily/day/{}".format(day)] = days[day]

    ordered_recent = [value for _, value in sorted(recent)] + legacy
    unique_recent = []
    for value in ordered_recent:
        if value not in unique_recent:
            unique_recent.append(value)
    for index, value in enumerate(unique_recent[:RECENT_CACHE_LIMIT]):
        desired["refs/ci/recent/{}".format(index)] = value

    _update_cache_refs(cache, desired, [ref for ref, _ in retained if ref not in desired])
    _command(("git", "--git-dir", str(cache), "gc", "--auto"))


def _split_bundle(bundle: Path, destination: Path) -> Tuple[List[Path], str]:
    size = bundle.stat().st_size
    if size < BUNDLE_PARTS:
        raise ValueError("source bundle is too small")
    digest = hashlib.sha256()
    parts = []
    with bundle.open("rb") as source:
        for index in range(BUNDLE_PARTS):
            amount = size * (index + 1) // BUNDLE_PARTS - size * index // BUNDLE_PARTS
            data = source.read(amount)
            digest.update(data)
            part = destination / "source.bundle.part.{:02d}".format(index)
            part.write_bytes(data)
            parts.append(part)
        if source.read(1):
            raise RuntimeError("unable to split source bundle")
    return parts, digest.hexdigest()


def _assemble_bundle(root: Path, expected: str) -> Path:
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("invalid bundle checksum")
    parts = [root / "source.bundle.part.{:02d}".format(index) for index in range(BUNDLE_PARTS)]
    if set(root.glob("source.bundle.part.*")) != set(parts) or not all(path.is_file() for path in parts):
        raise ValueError("incomplete source bundle")
    temporary = root / ".source.bundle.tmp"
    digest = hashlib.sha256()
    with temporary.open("wb") as output:
        for part in parts:
            with part.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
                    output.write(block)
    if digest.hexdigest() != expected:
        temporary.unlink()
        raise ValueError("source bundle checksum mismatch")
    bundle = root / "source.bundle"
    os.replace(str(temporary), str(bundle))
    for part in parts:
        part.unlink()
    return bundle


def remote_prepare(project: Path, run: Path) -> None:
    project = project.resolve()
    run = run.resolve()
    roots = tuple(Path(root).resolve() for root in load_config().remote.allowed_project_roots)
    if not _below_project_root(project, roots):
        raise ValueError("project must be below a configured project root")
    _run_relative(project, run)
    run_paths = tuple(run / name for name in ("source", "results", "jobs", "input"))
    if any(path.exists() for path in run_paths):
        raise ValueError("run directory already exists")
    for path in (*run_paths, project / "archives"):
        path.mkdir(parents=True, exist_ok=True)
    cache = _cache(project)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (cache.parent / "lock").open("w") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        if not cache.exists():
            _command(("git", "init", "--bare", str(cache)))
        _command(("git", "--git-dir", str(cache), "config", "fetch.unpackLimit", "1"))
        _cleanup_reservations(cache)
        retained = _retained_refs(cache)
        cached = _cached_commits(cache, retained)
        _reserve_cache(cache, run, retained)
        cleanup_archives(project)
    print(json.dumps({"cache_shas": sorted(set(cached))}))


def remote_receive(project: Path, run: Path, source_sha: str, bundle_checksum: str) -> None:
    if not SHA.fullmatch(source_sha):
        raise ValueError("invalid source SHA")
    cache = _cache(project.resolve())
    input_root = run.resolve() / "input"
    bundle = input_root / "source.bundle"
    if bundle_checksum != "-":
        bundle = _assemble_bundle(input_root, bundle_checksum)
    lock = cache.parent / "lock"
    with lock.open("w") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        if bundle_checksum != "-":
            heads = _command(("git", "bundle", "list-heads", str(bundle))).stdout.split()
            if len(heads) != 2 or heads != [source_sha, "HEAD"]:
                raise ValueError("bundle does not advertise the requested SHA")
            _command(("git", "--git-dir", str(cache), "bundle", "verify", str(bundle)))
            _command((
                "git", "--git-dir", str(cache), "fetch", str(bundle),
                "HEAD:refs/ci/{}".format(source_sha),
            ))
        else:
            if any(input_root.iterdir()):
                raise ValueError("unexpected bundle data for cache hit")
            _command(("git", "--git-dir", str(cache), "cat-file", "-e", source_sha + "^{commit}"))
        relative = _run_relative(project, run)
        _rotate_cache_refs(cache, "daily" if relative.parts[0] == "daily" else "recent", source_sha)
        source = run / "source"
        shutil.rmtree(str(source))
        source.mkdir()
        archive = subprocess.Popen(("git", "--git-dir", str(cache), "archive", source_sha), stdout=subprocess.PIPE)
        assert archive.stdout is not None
        extract = subprocess.run(("tar", "-xf", "-", "-C", str(source)), stdin=archive.stdout)
        archive.stdout.close()
        if archive.wait() or extract.returncode:
            raise RuntimeError("unable to extract source tree")
        _release_cache(cache, run)


def collect(run: Path) -> None:
    with tarfile.open(fileobj=sys.stdout.buffer, mode="w|gz") as archive:
        archive.add(str(run / "results"), arcname="results")


def cleanup_archives(project: Path) -> int:
    removed = 0
    root = project / "archives"
    now = time.time()
    if root.is_dir():
        for path in root.glob("*/*.tar.gz"):
            if now - path.stat().st_mtime > 72 * 3600:
                path.unlink()
                removed += 1
    return removed


def archive_run(project: Path, run: Path) -> Path:
    relative = run.resolve().relative_to((project / "runs").resolve())
    destination = project / "archives" / relative.parent / (relative.name + ".tar.gz")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp")
    with tarfile.open(str(temporary), "w:gz") as archive:
        archive.add(str(run / "results"), arcname="results")
        archive.add(str(run / "jobs"), arcname="jobs")
    os.replace(str(temporary), str(destination))
    shutil.rmtree(str(run))
    return destination


def _ssh(config: Path, target: str, command: Sequence[str], check: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(("ssh", "-F", str(config), target, shlex.join(command)), text=True, capture_output=True)
    if check and result.returncode:
        raise RuntimeError(result.stderr.strip() or "SSH command failed")
    return result


def _extract_results(stream: Any, destination: Path) -> None:
    members = 0
    size = 0
    with tarfile.open(fileobj=stream, mode="r|gz") as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            if path.is_absolute() or not path.parts or path.parts[0] != "results" or \
                    ".." in path.parts or not (member.isdir() or member.isfile()):
                raise ValueError("unsafe result archive member: {}".format(member.name))
            members += 1
            size += member.size
            if members > MAX_RESULT_MEMBERS or size > MAX_RESULT_BYTES:
                raise ValueError("result archive is too large")
            archive.extract(member, str(destination))


def _bundle_revision(repository: Path, cached: Sequence[str], source_sha: str) -> Optional[str]:
    if source_sha in cached:
        return None
    cached_set = set(cached)
    history = _command(("git", "rev-list", "--topo-order", "HEAD"), repository).stdout.splitlines()
    for base in history:
        if base in cached_set:
            return base + "..HEAD"
    return "HEAD"


def _validate_result(result: Any) -> Mapping[str, Any]:
    config = load_config()
    if type(result) is not dict or set(result) != {
        "protocol", "total", "passed", "failed", "infrastructure",
        "components", "cases",
    } or type(result["protocol"]) is not int or result["protocol"] != 1:
        raise ValueError("invalid result protocol")

    count_names = ("passed", "failed", "infrastructure", "total")
    counts = {name: result[name] for name in count_names}
    if any(type(value) is not int or value < 0 for value in counts.values()):
        raise ValueError("invalid result counts")

    components = result["components"]
    expected_components = [("build", "Compile")] + [
        (name, profile.label) for name, profile in config.resources.items()
    ]
    if type(components) is not list or len(components) != len(expected_components):
        raise ValueError("invalid result components")
    for item, identity in zip(components, expected_components):
        if type(item) is not dict or set(item) != {
            "name", "label", "state", "job_id", "slurm_state", "exit_code",
        } or any(type(item[name]) is not str for name in item) or \
                (item["name"], item["label"]) != identity:
            raise ValueError("invalid result component")

    rows = result["cases"]
    expected_cases = [
        case for resource in config.resources
        for case in config.cases if case.resource == resource
    ]
    if type(rows) is not list or len(rows) != len(expected_cases):
        raise ValueError("invalid result cases")
    for row, case in zip(rows, expected_cases):
        if type(row) is not dict or set(row) != {
            "case_id", "resource", "runner", "state", "exit_code",
            "slurm_state", "slurm_exit_code", "job_id", "elapsed_seconds",
        } or any(type(row[name]) is not str for name in (
            "case_id", "resource", "runner", "state", "slurm_state",
            "slurm_exit_code", "job_id",
        )) or (type(row["exit_code"]) is not int and row["exit_code"] is not None) or \
                type(row["elapsed_seconds"]) is not int or row["elapsed_seconds"] < 0 or \
                (row["case_id"], row["resource"], row["runner"]) != (
                    case.case_id, case.resource, case.runner,
                ) or row["state"] not in ("PASS", "FAIL", "TIMEOUT", "INFRA"):
            raise ValueError("invalid result case")

    passed = sum(row["state"] == "PASS" for row in rows)
    failed = sum(row["state"] in ("FAIL", "TIMEOUT") for row in rows)
    actual_counts = (passed, failed, len(rows) - passed - failed, len(rows))
    if tuple(counts[name] for name in count_names) != actual_counts:
        raise ValueError("inconsistent result counts")

    build_state = components[0]["state"]
    if build_state == "PASS":
        for component, resource in zip(components[1:], config.resources):
            states = [row["state"] for row in rows if row["resource"] == resource]
            expected = "PASS" if all(state == "PASS" for state in states) else (
                "FAIL" if any(state in ("FAIL", "TIMEOUT") for state in states) else "INFRA"
            )
            if component["state"] != expected:
                raise ValueError("inconsistent component state")
    elif build_state != "FAIL" or any(
        component["state"] != "SKIPPED" for component in components[1:]
    ) or any(row["state"] != "INFRA" for row in rows):
        raise ValueError("invalid build result")
    return result


def _read_result(path: Path) -> Mapping[str, Any]:
    if path.stat().st_size > MAX_REPORT_BYTES:
        raise ValueError("result report is too large")
    return _validate_result(json.loads(path.read_text(encoding="utf-8")))


def _print_progress(
    jobs: Sequence[Tuple[str, str]], states: Mapping[str, Mapping[str, int]],
    previous: Dict[str, Tuple[int, int, int]],
) -> None:
    for label, job in jobs:
        state = states[job]
        current = state["finished"], state["running"], state["total"]
        if previous.get(job) == current:
            continue
        previous[job] = current
        queued = state["total"] - state["finished"] - state["running"]
        details = ["{}/{} finished".format(state["finished"], state["total"])]
        if state["running"]:
            details.append("{} running".format(state["running"]))
        if queued:
            details.append("{} queued".format(queued))
        print("  {:<24} [{}] {}".format(label, job, ", ".join(details)), flush=True)


def _print_result(
    result: Optional[Mapping[str, Any]], artifacts: Path, remote_archive: str,
) -> None:
    root = artifacts.resolve()
    if result is None:
        print("\nGPU validation: INFRASTRUCTURE ERROR")
        print("No valid result report was produced.")
    else:
        state = "PASS" if result["passed"] == result["total"] else "FAIL"
        print("\nGPU validation: {}".format(state))
        print("{} passed, {} failed, {} infrastructure\n".format(
            result["passed"], result["failed"], result["infrastructure"],
        ))
        for component in _folder_components(result):
            print("  {:<24} {}".format(component["label"], component["state"]))
        print("\nSummary: {}".format(root / "results" / "summary.md"))
    print("Raw results: {}".format(root / "results"))
    print("Remote archive: {}".format(remote_archive))


def _upload_parts(parts: Sequence[Path], args: argparse.Namespace, run: Path) -> None:
    def upload(part: Path) -> None:
        _retry((
            "rsync", "-a", "--partial", "--info=stats2", "-e",
            "ssh -F {}".format(args.ssh_config), str(part),
            "{}:{}/input/".format(args.target, run),
        ))

    with ThreadPoolExecutor(max_workers=len(parts)) as pool:
        futures = [pool.submit(upload, part) for part in parts]
        for count, future in enumerate(as_completed(futures), 1):
            future.result()
            print("  Source upload: {}/{} parts transferred".format(count, len(parts)), flush=True)


def _artifact_path(args: argparse.Namespace) -> Path:
    if hasattr(args, "artifacts"):
        return args.artifacts
    return Path("/tmp") / "abacus_gpu_ci_{}".format(os.getuid()) / \
        args.namespace / "{}_{}".format(args.run_id, args.run_attempt)


def run(args: argparse.Namespace) -> int:
    args.ssh_config = args.ssh_config.expanduser()
    args.source_repository = args.source_repository.expanduser()
    automatic_artifacts = not hasattr(args, "artifacts")
    args.artifacts = _artifact_path(args).expanduser()
    repository = args.source_repository.resolve()
    config = load_config()
    if not all(NAME.fullmatch(value) for value in (args.target, args.namespace, args.run_id, args.run_attempt)):
        raise ValueError("invalid target or run name")
    actual = _command(("git", "rev-parse", "HEAD"), repository).stdout.strip()
    if args.source_sha == "HEAD":
        args.source_sha = actual
    if args.source_sha != actual or not SHA.fullmatch(actual):
        raise ValueError("source SHA must be the checked-out HEAD")
    requested = (args.project_root, *(str(root) for root in config.remote.allowed_project_roots))
    script = (
        "from pathlib import Path; import json,sys; "
        "print(json.dumps([str(Path(value).expanduser().resolve()) for value in sys.argv[1:]]))"
    )
    command = shlex.join(("python3", "-c", script, *requested))
    resolved = json.loads(_retry(("ssh", "-F", str(args.ssh_config), args.target, command)).stdout)
    if type(resolved) is not list or len(resolved) != len(requested) or \
            any(type(value) is not str for value in resolved):
        raise ValueError("invalid resolved project paths")
    project = Path(resolved[0])
    roots = tuple(Path(value) for value in resolved[1:])
    if not project.is_absolute() or any(not NAME.fullmatch(part) for part in project.parts[1:]):
        raise ValueError("project root must be a simple absolute path")
    if not _below_project_root(project, roots):
        raise ValueError("project root must be below a configured project root")
    args.project_root = str(project)
    run = project / "runs" / args.namespace / "{}-{}".format(args.run_id, args.run_attempt)
    print("Preparing remote run...", flush=True)
    if automatic_artifacts:
        artifact_root = args.artifacts.parents[1]
        artifact_root.mkdir(mode=0o700, exist_ok=True)
        artifact_root.chmod(0o700)
    args.artifacts.mkdir(mode=0o700, parents=True, exist_ok=True)
    _atomic(args.artifacts / "run.json", {
        "project_root": args.project_root, "run_root": str(run),
        "source_sha": args.source_sha,
    })
    remote_control = str(run / "control")
    _retry((
        "ssh", "-F", str(args.ssh_config), args.target,
        shlex.join(("mkdir", "-p", remote_control)),
    ))
    _retry((
        "rsync", "-a", "--delete", "--info=stats2", "-e",
        "ssh -F {}".format(args.ssh_config), str(ROOT) + "/",
        "{}:{}/".format(args.target, remote_control),
    ))
    prepared = _ssh(args.ssh_config, args.target, (
        "python3", remote_control + "/runner.py", "prepare", args.project_root, str(run),
    ))
    data = json.loads(prepared.stdout)
    if type(data) is not dict or set(data) != {"cache_shas"} or type(data["cache_shas"]) is not list or \
            any(type(value) is not str or not SHA.fullmatch(value) for value in data["cache_shas"]) or \
            len(set(data["cache_shas"])) != len(data["cache_shas"]):
        raise ValueError("invalid remote source cache")
    print("Uploading source...", flush=True)
    with tempfile.TemporaryDirectory() as directory:
        bundle = Path(directory) / "source.bundle"
        revision = _bundle_revision(repository, data["cache_shas"], args.source_sha)
        checksum = "-"
        if revision:
            _command(("git", "bundle", "create", str(bundle), revision), repository)
            parts, checksum = _split_bundle(bundle, Path(directory))
            _upload_parts(parts, args, run)
        else:
            print("  Source already cached", flush=True)
    _ssh(args.ssh_config, args.target, (
        "python3", remote_control + "/runner.py", "receive",
        args.project_root, str(run), args.source_sha, checksum,
    ))
    print("Building and testing...", flush=True)
    launch = "nohup python3 {} remote-run {} > {}/results/coordinator.log 2>&1 < /dev/null &".format(
        shlex.quote(remote_control + "/runner.py"), shlex.quote(str(run)), shlex.quote(str(run))
    )
    _ssh(args.ssh_config, args.target, ("bash", "-lc", launch))
    failures = 0
    printed_lines = 0
    poll_seconds = config.poll_seconds
    marker = "\n---DONE---\n"
    while True:
        log_path = shlex.quote(str(run / "results" / "coordinator.log"))
        done_path = shlex.quote(str(run / "results" / "done.json"))
        read_status = "cat {0} 2>/dev/null || true; printf '\\n---DONE---\\n'; cat {1} 2>/dev/null || true".format(
            done_path, log_path,
        )
        status = _ssh(args.ssh_config, args.target, (
            "bash", "-lc", read_status,
        ), check=False)
        done_text, _, progress_text = status.stdout.partition(marker)
        if status.returncode:
            failures += 1
            if failures == 10:
                raise RuntimeError("lost contact with the remote cluster")
        else:
            failures = 0
            if progress_text.strip():
                lines = progress_text.splitlines()
                for line in lines[printed_lines:]:
                    print(line, flush=True)
                printed_lines = len(lines)
        if not status.returncode and done_text.strip():
            done = json.loads(done_text)
            if type(done) is not dict or set(done) != {"returncode"} or \
                    type(done["returncode"]) is not int or done["returncode"] not in (0, 1, 2):
                raise ValueError("invalid remote completion status")
            break
        time.sleep(poll_seconds)
    print("Downloading results...", flush=True)
    command = (
        "ssh", "-F", str(args.ssh_config), args.target,
        shlex.join(("python3", remote_control + "/runner.py", "collect", str(run))),
    )
    archive = args.artifacts / ".results.tar.gz"
    _retry_download(command, archive)
    with archive.open("rb") as stream:
        _extract_results(stream, args.artifacts)
    archive.unlink()
    result_path = args.artifacts / "results" / "result.json"
    result: Optional[Mapping[str, Any]] = None
    if result_path.is_file():
        result = _read_result(result_path)
        expected_returncode = 0 if result["passed"] == result["total"] else 1
    else:
        expected_returncode = 2
    if done["returncode"] != expected_returncode:
        raise ValueError("completion status does not match result")
    print("Archiving remote run...", flush=True)
    archived = project / "archives" / args.namespace / (run.name + ".tar.gz")
    archive_command = "test -f {archive} || python3 {runner} archive {project} {run} >/dev/null; test -f {archive}".format(
        archive=shlex.quote(str(archived)), runner=shlex.quote(remote_control + "/runner.py"),
        project=shlex.quote(str(project)), run=shlex.quote(str(run)),
    )
    _retry(("ssh", "-F", str(args.ssh_config), args.target, archive_command))
    _print_result(result, args.artifacts, str(archived))
    return done["returncode"]


def report(args: argparse.Namespace) -> int:
    if not args.result.is_file():
        components = [{"name": "infrastructure", "label": "Infrastructure", "state": "INFRA"}]
        values = {"available": "false", "passed": "", "failed": "", "infrastructure": "", "total": ""}
    else:
        result = _read_result(args.result)
        components = [
            {key: item[key] for key in ("name", "label", "state")}
            for item in _folder_components(result)
        ]
        counts = {name: result[name] for name in ("passed", "failed", "infrastructure", "total")}
        values = {"available": "true", **{name: str(value) for name, value in counts.items()}}
        if args.summary:
            args.summary.write_text(_result_markdown(result), encoding="utf-8")
    values["components"] = json.dumps(components, separators=(",", ":"))
    with args.output.open("a", encoding="utf-8") as stream:
        for name, value in values.items():
            stream.write("{}={}\n".format(name, value))
    return 0


def _gh(path: str, method: str = "GET", fields: Optional[Mapping[str, str]] = None) -> Any:
    command = ["gh", "api", "--method", method, path]
    for name, value in (fields or {}).items():
        command.extend(("-f", "{}={}".format(name, value)))
    return json.loads(_command(command).stdout)


def github_admit() -> int:
    event = os.environ["GITHUB_EVENT_NAME"]
    repository = os.environ["GITHUB_REPOSITORY"]
    values = {"accepted": "true", "pr_number": "", "check_id": "", "comment_id": ""}
    if event == "schedule":
        upstream = os.environ.get("UPSTREAM_REPOSITORY", "deepmodeling/abacus-develop")
        metadata = _gh("repos/{}".format(upstream))
        commit = _gh("repos/{}/commits/{}".format(upstream, metadata["default_branch"]))
        values.update(source_repository=upstream, source_sha=commit["sha"], namespace="daily")
    elif event == "workflow_dispatch":
        values.update(
            source_repository=repository,
            source_sha=os.environ["MANUAL_SOURCE_SHA"], namespace="manual",
        )
    else:
        event_data = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8"))
        user = event_data["comment"]["user"]["login"]
        permission = _gh("repos/{}/collaborators/{}/permission".format(repository, user))
        if permission.get("role_name") not in ("admin", "maintain", "write", "triage"):
            raise ValueError("commenter needs Triage permission")
        number = str(event_data["issue"]["number"])
        pull = _gh("repos/{}/pulls/{}".format(repository, number))
        if pull["state"] != "open":
            raise ValueError("pull request is not open")
        values.update(
            source_repository=pull["head"]["repo"]["full_name"],
            source_sha=pull["head"]["sha"], namespace="pr-" + number,
            pr_number=number,
        )
        if not SHA.fullmatch(values["source_sha"]):
            raise ValueError("invalid source SHA")
        comment = _gh("repos/{}/issues/{}/comments".format(repository, number), "POST", {
            "body": (
                "## GPU validation: queued\n\n"
                "[Open the Actions run]({})\n\n"
                "Candidate: `{}`"
            ).format(os.environ["RUN_URL"], values["source_sha"]),
        })
        values["comment_id"] = str(comment["id"])
        check = _gh("repos/{}/check-runs".format(repository), "POST", {
            "name": "GPU validation", "head_sha": pull["head"]["sha"],
            "status": "in_progress",
        })
        values["check_id"] = str(check["id"])
    if not SHA.fullmatch(values["source_sha"]):
        raise ValueError("invalid source SHA")
    with Path(os.environ["GITHUB_OUTPUT"]).open("a", encoding="utf-8") as stream:
        for name, value in values.items():
            stream.write("{}={}\n".format(name, value))
    return 0


def github_finish() -> int:
    repository = os.environ["GITHUB_REPOSITORY"]
    result = os.environ["GPU_RESULT"]
    conclusion = "success" if result == "success" else "failure"
    errors = []
    check_id = os.environ["CHECK_ID"]
    if check_id:
        try:
            _gh("repos/{}/check-runs/{}".format(repository, check_id), "PATCH", {
                "status": "completed", "conclusion": conclusion,
            })
        except (OSError, RuntimeError, ValueError) as error:
            errors.append(error)
    pr = os.environ["PR_NUMBER"]
    if pr:
        body = (
            "## GPU validation: {}\n\n"
            "GPU cases: **{} passed, {} failed, {} infrastructure**.\n\n"
            "[Open the Actions run]({}) | [Download raw test files]({})\n\n"
            "Candidate: `{}`"
        ).format(
            result, os.environ.get("GPU_PASSED") or "?", os.environ.get("GPU_FAILED") or "?",
            os.environ.get("GPU_INFRASTRUCTURE") or "?", os.environ["RUN_URL"],
            os.environ.get("ARTIFACT_URL") or os.environ["RUN_URL"], os.environ["SOURCE_SHA"],
        )
        try:
            _gh("repos/{}/issues/comments/{}".format(repository, os.environ["COMMENT_ID"]), "PATCH", {"body": body})
        except (OSError, RuntimeError, ValueError) as error:
            errors.append(error)
    if errors:
        raise RuntimeError("; ".join(str(error) for error in errors))
    return 0


def parser() -> argparse.ArgumentParser:
    main = argparse.ArgumentParser(description=__doc__)
    commands = main.add_subparsers(
        title="commands", dest="command", required=True, metavar="{run,config}",
    )
    client = commands.add_parser(
        "run", help="build and run the GPU matrix on a remote cluster",
        description="Transfer one committed ABACUS revision, build it, and run the GPU matrix through Slurm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    client.add_argument(
        "--ssh-config", type=Path, default=Path("~/.ssh/config"),
        help="SSH config file containing the target host alias",
    )
    client.add_argument(
        "--target", default="gpu-ci",
        help="host alias in the SSH config",
    )
    client.add_argument(
        "--project-root", default=load_config().remote.project_root,
        help="remote directory below a configured project root for caches, runs, and archives",
    )
    client.add_argument(
        "--source-repository", type=Path, default=REPOSITORY_ROOT,
        help="local ABACUS Git checkout to transfer",
    )
    client.add_argument(
        "--source-sha", default="HEAD",
        help="commit to test; it must resolve to the checkout's HEAD",
    )
    client.add_argument(
        "--namespace", default="manual",
        help="group for the remote run, such as manual, daily, or pr-123",
    )
    client.add_argument(
        "--run-id", default=str(int(time.time())),
        help="unique identifier for this run",
    )
    client.add_argument(
        "--run-attempt", default="1",
        help="attempt number within the run ID",
    )
    client.add_argument(
        "--artifacts", type=Path, default=argparse.SUPPRESS,
        help=(
            "local directory for downloaded results and client logs "
            "(default root: /tmp/abacus_gpu_ci_<uid>; "
            "run directory: <namespace>/<run_id>_<attempt>)"
        ),
    )
    commands.add_parser(
        "config", help="validate and summarize config.ini",
        description="Validate config.ini and print its remote, resource, and case configuration.",
    )

    def internal(name: str) -> argparse.ArgumentParser:
        return commands.add_parser(name)

    remote = internal("remote-run")
    remote.add_argument("run", type=Path)
    worker_parser = internal("worker")
    for name in ("source", "control", "install", "results", "manifest"):
        worker_parser.add_argument(name, type=Path)
    prepare = internal("prepare")
    prepare.add_argument("project", type=Path)
    prepare.add_argument("run", type=Path)
    receive = internal("receive")
    receive.add_argument("project", type=Path)
    receive.add_argument("run", type=Path)
    receive.add_argument("source_sha")
    receive.add_argument("bundle_checksum")
    collect_parser = internal("collect")
    collect_parser.add_argument("run", type=Path)
    archive = internal("archive")
    archive.add_argument("project", type=Path)
    archive.add_argument("run", type=Path)
    report_parser = internal("report")
    report_parser.add_argument("--result", required=True, type=Path)
    report_parser.add_argument("--output", required=True, type=Path)
    report_parser.add_argument("--summary", type=Path)
    internal("github-admit")
    internal("github-finish")
    return main


def main() -> int:
    if sys.version_info < (3, 8):
        raise RuntimeError("Python 3.8 or newer is required")
    args = parser().parse_args()
    if args.command == "config":
        config = load_config()
        print(json.dumps({
            "site": {
                "name": config.site.name, "url": config.site.url,
                "acknowledgement": config.site.acknowledgement,
            },
            "remote": {
                "host": config.remote.host, "port": config.remote.port,
                "user": config.remote.user, "project_root": config.remote.project_root,
                "allowed_project_roots": list(map(str, config.remote.allowed_project_roots)),
            },
            "resources": list(config.resources), "cases": len(config.cases),
        }))
        return 0
    if args.command == "remote-run":
        return remote_run(args.run)
    if args.command == "worker":
        return worker(args.source, args.control, args.install, args.results, args.manifest)
    if args.command == "prepare":
        remote_prepare(args.project, args.run)
    elif args.command == "receive":
        remote_receive(args.project, args.run, args.source_sha, args.bundle_checksum)
    elif args.command == "collect":
        collect(args.run)
    elif args.command == "archive":
        print(archive_run(args.project, args.run))
    elif args.command == "run":
        return run(args)
    elif args.command == "report":
        return report(args)
    elif args.command == "github-admit":
        return github_admit()
    elif args.command == "github-finish":
        return github_finish()
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError, SlurmError) as error:
        print("gpu-ci: {}".format(error), file=sys.stderr)
        sys.exit(2)
