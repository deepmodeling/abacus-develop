"""Small Slurm adapter for the GPU matrix."""

import re
import subprocess
import time
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple


Terminal = Tuple[str, str]
TERMINAL = {
    "BOOT_FAIL", "CANCELLED", "COMPLETED", "DEADLINE", "FAILED",
    "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT",
}


class SlurmError(RuntimeError):
    pass


class Slurm:
    def __init__(self, poll_seconds: int = 10) -> None:
        self.poll_seconds = poll_seconds
        self.jobs: Dict[str, Optional[int]] = {}

    @staticmethod
    def _run(command: Sequence[str]) -> str:
        result = subprocess.run(command, text=True, capture_output=True)
        if result.returncode:
            raise SlurmError(result.stderr.strip() or result.stdout.strip())
        return result.stdout

    def submit(self, script: Path, array_count: Optional[int] = None) -> str:
        output = self._run(("sbatch", "--parsable", str(script))).strip()
        match = re.fullmatch(r"([0-9]+)(?:;[A-Za-z0-9_.-]+)?", output)
        if not match:
            raise SlurmError("invalid sbatch output: {!r}".format(output))
        job = match.group(1)
        self.jobs[job] = array_count
        return job

    def wait(self, jobs: Sequence[str]) -> Dict[str, Terminal]:
        ids = ",".join(jobs)
        failures = 0
        while True:
            try:
                if not self._run(("squeue", "--noheader", "--jobs=" + ids)).strip():
                    break
                failures = 0
            except SlurmError:
                failures += 1
                if failures == 6:
                    raise
            time.sleep(self.poll_seconds)

        for _ in range(30):
            rows = self._accounting(ids)
            required = []
            for job in jobs:
                count = self.jobs[job]
                required.extend(
                    [job] if count is None else
                    ["{}_{}".format(job, index) for index in range(count)]
                )
            if all(job in rows for job in required):
                return {job: rows[job] for job in required}
            time.sleep(self.poll_seconds)
        raise SlurmError("Slurm accounting did not become complete")

    def _accounting(self, jobs: str) -> Dict[str, Terminal]:
        output = self._run((
            "sacct", "--noheader", "--allocations", "--parsable2",
            "--jobs=" + jobs, "--format=JobID,State,ExitCode",
        ))
        rows: Dict[str, Terminal] = {}
        for line in output.splitlines():
            fields = line.strip().split("|")
            if len(fields) != 3 or "." in fields[0]:
                continue
            state = fields[1].split()[0].rstrip("+") if fields[1].strip() else ""
            if re.fullmatch(r"[0-9]+(?:_[0-9]+)?", fields[0]) and state in TERMINAL:
                rows[fields[0]] = state, fields[2]
        return rows

    def cancel(self) -> None:
        if self.jobs:
            subprocess.run(("scancel", *self.jobs), check=False)
