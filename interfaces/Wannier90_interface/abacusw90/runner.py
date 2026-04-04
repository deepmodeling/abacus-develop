"""
Runner utility for executing external processes (ABACUS, Wannier90).
"""

import subprocess
import os
from pathlib import Path
from typing import Optional

def run_command(
    cmd: str, 
    cwd: str, 
    log_file: Optional[str] = None,
    env: Optional[dict] = None,
    shell: bool = True
) -> int:
    """
    Executes a shell command in a specific directory.
    
    Args:
        cmd (str): Command string to execute.
        cwd (str): Working directory.
        log_file (str, optional): File to write stdout/stderr. Defaults to None.
        env (dict, optional): Environment variables. Defaults to None.
        shell (bool, optional): Use shell. Defaults to True.
        
    Returns:
        int: Return code of the process.
    """
    print(f">>> Executing: {cmd} (in {cwd})")
    
    # Prepare environment
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
        
    # Open log file if specified
    log_f = open(log_file, 'w') if log_file else None
    
    try:
        process = subprocess.Popen(
            cmd,
            cwd=cwd,
            shell=shell,
            env=run_env,
            stdout=log_f if log_f else subprocess.PIPE,
            stderr=log_f if log_f else subprocess.PIPE,
            text=True
        )
        
        process.wait()
        
        if process.returncode != 0:
            print(f"Warning: Command '{cmd}' exited with return code {process.returncode}")
            
        return process.returncode
        
    except Exception as e:
        print(f"Error executing command: {e}")
        return 1
        
    finally:
        if log_f:
            log_f.close()

def check_file_exists(filepath: Path, raise_error: bool = True) -> bool:
    """Check if a file exists."""
    if not filepath.exists():
        if raise_error:
            raise FileNotFoundError(f"Required file not found: {filepath}")
        return False
    return True
