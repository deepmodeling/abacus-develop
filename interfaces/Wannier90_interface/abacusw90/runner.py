"""runner.py - Runner utility for executing external processes (ABACUS, Wannier90)."""

import subprocess
import os
from pathlib import Path
from typing import Optional


def run_command(
    cmd: str,
    cwd: str,
    log_file: Optional[str] = None,
    env: Optional[dict] = None,
    shell: bool = True,
    label: str = "command",
    check: bool = True,
) -> int:
    """Executes a shell command in a specific directory."""
    print(f">>> [{label}] Executing: {cmd}")
    print(f">>> [{label}] Working directory: {cwd}")

    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    if log_file is None:
        log_file = str(Path(cwd) / f".{label}.log")

    try:
        with open(log_file, "w") as log_f:
            process = subprocess.Popen(
                cmd,
                cwd=cwd,
                shell=shell,
                env=run_env,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
            )
            process.wait()

        if process.returncode != 0:
            with open(log_file, "r") as f:
                lines = f.readlines()
            tail = "".join(lines[-20:]) if lines else "(empty)"
            msg = (
                f"[{label}] Command failed with return code {process.returncode}.\n"
                f"  Command : {cmd}\n"
                f"  CWD     : {cwd}\n"
                f"  Log file: {log_file}\n"
                f"  Last 20 lines:\n{tail}"
            )
            if check:
                raise RuntimeError(msg)
            else:
                print(f"Warning: {msg}")
        else:
            print(f">>> [{label}] Finished successfully (rc=0)")

        return process.returncode

    except FileNotFoundError:
        raise RuntimeError(
            f"[{label}] Executable not found in command: '{cmd}'.\n"
            f"  Check that '{cmd.split()[0]}' is installed and in PATH."
        )
    except Exception as e:
        raise RuntimeError(f"[{label}] Error executing command: {e}") from e


def check_file_exists(
    filepath,
    raise_error: bool = True,
    hint: str = "",
) -> bool:
    """Check if a file exists, with diagnostic context on failure."""
    filepath = Path(filepath)
    if not filepath.exists():
        if raise_error:
            parent = filepath.parent
            if parent.exists():
                siblings = [f.name for f in parent.iterdir()]
                dir_info = f"Directory contents: {siblings}"
            else:
                dir_info = f"Parent directory does not exist: {parent}"

            msg = f"Required file not found: {filepath}"
            if hint:
                msg += f"\n  Hint: {hint}"
            msg += f"\n  {dir_info}"
            raise FileNotFoundError(msg)
        return False
    return True
