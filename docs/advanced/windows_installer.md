# Windows One-Click Installer (WSL2 + conda-forge)

A lightweight installer that brings ABACUS to Windows via WSL2 and
conda-forge. No C++ toolchain, no MPI build, no manual dependency juggling —
run `install-abacus.bat` once and type `abacus` from any Windows terminal.

The scripts live in the repository under `tools/windows/`. This page is the
user-facing documentation for the same scripts.

## How it works

ABACUS depends on a heavy Linux-native scientific stack (OpenMPI, ScaLAPACK,
ELPA, FFTW, libxc, OpenBLAS, …) that is painful to build natively on Windows.
Instead of porting, this installer provisions a standard Linux environment
inside WSL2 and exposes it through thin Windows launchers.

The pipeline, end to end:

1. **`install-abacus.bat`** (runs on Windows, requires admin)
   - Checks the Windows build (≥ 19041) and whether WSL is installed; if not,
     runs `wsl --install --no-launch` and asks the user to reboot once.
   - Optionally enables TUNA (Tsinghua) mirrors for users in Mainland China.
   - Prompts for an ABACUS version (blank = latest on conda-forge; an exact
     version like `3.7.4` is pinned; a match-spec like `>=3.7,<3.8` is passed
     through to conda).
   - Detects the target distribution (`Ubuntu-22.04`) by querying the WSL
     registry key `HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss`. This
     is immune to UTF-16 parsing pitfalls and to Store appx leftovers that
     can make `wsl -d <name> -- true` falsely report success.
   - Calls `wsl --install -d Ubuntu-22.04 --no-launch` if the distro is
     missing, then verifies the registry entry appeared.
   - Invokes `provision.sh` inside the distribution, stripping any `\r`
     bytes on the fly (`sed 's/\r$//' script | bash`) so Windows line
     endings don't break shell parsing.
   - Writes two small launcher `.cmd` files and adds them to the user PATH
     via PowerShell (avoiding `setx`'s 1024-character truncation).

2. **`provision.sh`** (runs as root inside the WSL distribution)
   - Optionally rewrites `/etc/apt/sources.list` to TUNA.
   - `apt-get install`s a minimal set of prerequisites (curl, ca-certificates,
     bzip2).
   - Downloads the Miniforge installer (from GitHub or TUNA's GitHub-releases
     mirror) and installs it to `/opt/abacus-miniforge`.
   - `conda create -n abacus_env -c conda-forge abacus` (or the TUNA
     conda-forge channel) — a single package pulls in the entire scientific
     runtime. conda-forge ships `abacus` for `linux-64` and `linux-aarch64`,
     which is exactly what WSL2 provides.
   - Writes two system-wide launchers, `/usr/local/bin/abacus` and
     `/usr/local/bin/abacus-mpi`, that activate the env and exec the real
     binary. Both set `OMP_NUM_THREADS=1` by default to avoid thread
     oversubscription. `abacus-mpi` additionally sets OpenMPI 4/5 + PRRTE
     "allow run as root" environment variables and passes
     `--allow-run-as-root` to `mpirun`, so the default WSL root user can
     launch parallel jobs without creating a non-root user.

3. **Windows launchers** (`abacus.cmd`, `abacus-mpi.cmd`)
   - Added to `%LOCALAPPDATA%\ABACUS\bin` and the user PATH.
   - Each launcher sets `WSLENV=OMP_NUM_THREADS:MKL_NUM_THREADS:OPENBLAS_NUM_THREADS:...`
     so thread-count overrides set on the Windows side are visible inside WSL.
   - Body is just:
     ```
     wsl -d Ubuntu-22.04 --cd "%CD%" -- abacus %*
     ```
     `--cd "%CD%"` maps the current Windows directory (`C:\…\case`) to its
     WSL path (`/mnt/c/…/case`), so users can `cd` into a case directory in
     `cmd`/PowerShell/Terminal and just type `abacus`.

4. **`uninstall-abacus.bat`**
   - Reads `install-state.txt` (written by the installer) to learn whether
     the Ubuntu-22.04 distribution was pre-existing or added by us.
   - If it was added by us, prompts whether to `wsl --unregister` the entire
     distribution, or to only wipe `/opt/abacus-miniforge` and the launchers.
   - If it was pre-existing, only the ABACUS files inside are removed.
   - Cleans Windows-side launchers and removes the bin directory from the
     user PATH.
   - Does **not** touch WSL itself (runtime, Windows optional features, or
     other distributions). See *Uninstallation* below for how to fully
     remove WSL if you want to.

## Requirements

- Windows 10 build 19041 (2004) or newer, or any Windows 11.
- Administrator privileges for the first run (to enable WSL features).
- Virtualization enabled in BIOS/UEFI.
- ~2 GB free disk space (Ubuntu + conda env).
- Network access to GitHub and conda-forge, or to TUNA if you choose the
  China mirror option.

## Installation

1. Clone or download this repository.
2. In `tools/windows/`, right-click `install-abacus.bat` → **Run as administrator**.
3. Answer the China-mirror prompt (`y` recommended inside Mainland China).
   Then pick an ABACUS version when prompted (leave blank for the latest on
   conda-forge; for a pinned install, type an exact version such as `3.7.4`).
4. If this is the first time WSL is installed on the machine, the script
   will ask you to reboot and run it again.
5. Wait for `[*] Provisioning ABACUS …` to finish (5–15 minutes on first
   run; most of it is the conda-forge download).
6. When you see `Installation complete!`, **open a new terminal window**
   (so the updated PATH takes effect) and verify:
   ```
   abacus --version
   ```

## Usage

Serial run in any case directory:

```
cd path\to\my_case
abacus
```

Parallel run with 4 MPI ranks:

```
abacus-mpi -n 4
```

Hybrid MPI + OpenMP (e.g. 4 MPI ranks × 2 threads each):

```
set OMP_NUM_THREADS=2
abacus-mpi -n 4
```

Set `OMP_NUM_THREADS` (and/or `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`) in
your Windows shell and the launcher will forward the value into WSL through
`WSLENV`. Unset, it defaults to 1 — a safe choice when running pure MPI.

Interactive Linux shell (for advanced debugging, manually running
`mpirun`, inspecting logs, etc.):

```
wsl -d Ubuntu-22.04
```

Inside the shell you can `conda activate abacus_env` to get access to
`mpirun`, `mpiexec`, and other tools from the conda environment.

## Uninstallation

### Standard: remove ABACUS only

Run `uninstall-abacus.bat`. This handles the common case:

- Removes `/opt/abacus-miniforge` and the `abacus` / `abacus-mpi` launchers
  inside WSL.
- If the installer added the `Ubuntu-22.04` distribution, asks whether you
  also want to `wsl --unregister` it (pick `y` to reclaim the disk space,
  `n` to keep the Linux environment for other uses).
- Deletes `%LOCALAPPDATA%\ABACUS\` and removes its `bin\` directory from
  your user PATH.

This is enough for almost every user. WSL itself and any *other* WSL
distributions you have stay untouched — important because WSL is commonly
shared with Docker Desktop, VS Code Remote, and other toolchains.

### Nuclear: remove WSL itself

Only do this if you truly have no other use for WSL on this machine.
Removing WSL will break Docker Desktop, VS Code Remote-WSL, any other
Linux distros you have, and so on. Run the following in an elevated
PowerShell:

```powershell
# 1. Unregister every WSL distribution (this wipes all their files).
wsl --list --quiet | ForEach-Object { wsl --unregister $_.Trim() }

# 2. Uninstall the WSL runtime itself (wsl.exe + the Linux kernel package).
#    This is Microsoft's official command; it does not disable the Windows
#    optional features and does not touch any distro appx packages.
wsl --uninstall

# 3. Optionally remove leftover distribution appx packages from the Store
#    (e.g. "Ubuntu 22.04 LTS"). `wsl --unregister` deletes the data only;
#    the Store app that installs the distro is separate.
Get-AppxPackage *Ubuntu* | Remove-AppxPackage

# 4. Optional: disable the Windows optional features (requires a reboot).
#    Skip this step if anything else on the machine still uses Hyper-V
#    virtualization (Docker Desktop, Windows Sandbox, Hyper-V VMs, ...).
dism.exe /online /disable-feature /featurename:Microsoft-Windows-Subsystem-Linux /norestart
dism.exe /online /disable-feature /featurename:VirtualMachinePlatform /norestart

# 5. Remove user config if present.
Remove-Item "$env:UserProfile\.wslconfig" -ErrorAction SilentlyContinue

# 6. Reboot to finalize.
Restart-Computer
```

After the reboot `wsl.exe` no longer exists. If you also ran step 4, the
Hyper-V virtualization layer used by WSL2 is disabled.

> On older Windows builds where `wsl --uninstall` is not available (WSL
> shipped via the in-box `wsl.exe` rather than the Store package), use
> `Get-AppxPackage *WindowsSubsystemForLinux* | Remove-AppxPackage` as a
> fallback for step 2.

## Performance notes

- Files under `/mnt/c/...` are served through the 9P protocol and are
  noticeably slower than native ext4. For heavy I/O (large SCF, MD
  trajectories), run the case from inside the WSL filesystem:
  ```
  wsl -d Ubuntu-22.04
  cp -r /mnt/c/path/to/case ~/case
  cd ~/case
  abacus
  ```
- The first `wsl` invocation after a boot triggers a 10–30 s VM cold start.
- OpenMPI runs all ranks inside a single WSL2 VM, so there is no network
  overhead between ranks — you get near-native parallel performance.

## File layout

```
tools/windows/
├── install-abacus.bat     # Windows entry point (admin, interactive)
├── uninstall-abacus.bat   # Clean removal, optionally including the distro
├── provision.sh           # Linux-side installer (runs as root in WSL)
├── .gitattributes         # Pin *.sh to LF, *.bat/*.cmd to CRLF
└── README.md              # Mirror of this page, shipped with the scripts
```

Artifacts created at install time:

```
%LOCALAPPDATA%\ABACUS\
├── bin\
│   ├── abacus.cmd         # Windows launcher (serial)
│   └── abacus-mpi.cmd     # Windows launcher (MPI)
└── install-state.txt      # Records whether we created the WSL distro

Inside WSL (Ubuntu-22.04):
/opt/abacus-miniforge/                       # Private Miniforge install
/opt/abacus-miniforge/envs/abacus_env/       # conda env holding abacus
/usr/local/bin/abacus, /usr/local/bin/abacus-mpi   # Linux-side launchers
```

## Design choices and trade-offs

- **Why WSL2 + conda-forge instead of a native Windows build?** ABACUS's
  MPI + ScaLAPACK + ELPA stack has no reliable native Windows build. Going
  through WSL2 lets us reuse the Linux binaries conda-forge already ships,
  turning a multi-week porting problem into a 200-line shell script.
- **Why a dedicated `Ubuntu-22.04` distribution?** conda-forge ABACUS is
  built against glibc from 22.04-era Ubuntu. Using `Ubuntu` (rolling) risks
  mismatches; pinning the version keeps the install reproducible.
- **Why put conda under `/opt/abacus-miniforge` rather than `/root`?**
  Clean uninstall path, clear ownership, and doesn't interfere with a
  user's personal conda install if they later add one inside the same
  distribution.
- **Why not ship a pre-built WSL rootfs?** Would cut first-run time from
  ~10 min to ~1 min, but balloons the installer from a few KB of scripts
  to 300–500 MB, requires CI infrastructure, and needs code signing to
  avoid SmartScreen warnings. A scripted online installer is the
  lowest-maintenance starting point; the pre-built rootfs path remains
  open for a future v1.
