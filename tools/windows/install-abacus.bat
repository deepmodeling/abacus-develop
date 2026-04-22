@echo off
setlocal EnableExtensions EnableDelayedExpansion
title ABACUS for Windows (WSL2) Installer

set "DISTRO=Ubuntu-22.04"
set "INSTALL_DIR=%LOCALAPPDATA%\ABACUS"
set "BIN_DIR=%INSTALL_DIR%\bin"

echo =========================================
echo   ABACUS Windows Installer (WSL2)
echo =========================================
echo.

REM --- 1. Admin check ---
net session >nul 2>&1
if errorlevel 1 (
    echo [!] Please right-click this script and choose "Run as administrator".
    pause
    exit /b 1
)

REM --- 2. Windows build check (WSL2 requires 19041+) ---
for /f "tokens=4-6 delims=. " %%a in ('ver') do set "BUILD=%%c"
if %BUILD% LSS 19041 (
    echo [!] Windows build %BUILD% detected.
    echo [!] WSL2 requires Windows 10 build 19041 or newer, or Windows 11.
    pause
    exit /b 1
)

REM --- 2b. Offer China mirror (TUNA) for faster downloads ---
if not defined ABACUS_CHINA_MIRROR (
    set "ABACUS_CHINA_MIRROR=0"
    set /p CHINA_INPUT="Use TUNA mirror for faster downloads (recommended inside Mainland China)? (y/N): "
    if /i "!CHINA_INPUT!"=="y" set "ABACUS_CHINA_MIRROR=1"
)

REM --- 3. Ensure WSL itself is present ---
where wsl >nul 2>&1
if errorlevel 1 (
    echo [*] WSL not found. Installing WSL2 runtime...
    wsl --install --no-launch
    if errorlevel 1 (
        echo [!] WSL installation failed.
        echo [!] Make sure virtualization is enabled in BIOS/UEFI.
        pause & exit /b 1
    )
    echo.
    echo [!] WSL was just installed. Please REBOOT Windows,
    echo [!] then run this script again to continue.
    pause
    exit /b 0
)

REM --- 4. Ensure target distro registered ---
echo [*] Updating WSL2 runtime if needed (Microsoft's progress bar shows below)...
wsl --update
wsl --set-default-version 2

REM Detect distro via registry (reliable, ANSI output, no UTF-16 pitfalls).
set "DISTRO_PREEXISTED=0"
reg query "HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss" /s /v DistributionName 2>nul | findstr /i /c:"%DISTRO%" >nul
if not errorlevel 1 set "DISTRO_PREEXISTED=1"

REM NOTE: keep this block flat (no nested parens) so `if errorlevel 1` after
REM each wsl.exe call reads the *runtime* errorlevel, not the parse-time one.
if "%DISTRO_PREEXISTED%"=="1" goto :distro_ready

echo [*] Installing %DISTRO% into WSL ^(Microsoft's download progress shows below^)...
wsl --install -d %DISTRO% --no-launch
if errorlevel 1 (
    echo [!] Failed to install %DISTRO%. Reboot and retry, or run
    echo [!]     wsl --install -d %DISTRO%
    echo [!] manually and finish any setup prompt, then re-run this script.
    pause
    exit /b 1
)

echo [*] Initializing %DISTRO% (first cold-start can take 30-60s)...
wsl -d %DISTRO% -u root -- /bin/true >nul 2>&1

reg query "HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss" /s /v DistributionName 2>nul | findstr /i /c:"%DISTRO%" >nul
if errorlevel 1 (
    echo [!] %DISTRO% still not registered after install. Possible causes:
    echo [!]   - Microsoft Store blocked/unreachable
    echo [!]   - A reboot is required to finish the initial WSL setup
    echo [!] Reboot, then re-run this script. Or run manually:
    echo [!]     wsl --install -d %DISTRO%
    pause
    exit /b 1
)

:distro_ready

REM --- 5. Run provisioning script inside the distro ---
echo [*] Provisioning ABACUS via conda-forge inside %DISTRO%.
echo     This downloads ~400 MB and takes 5-15 minutes on first run.
echo.

echo [*] Translating script path into WSL...
set "WSL_SCRIPT="
for /f "usebackq delims=" %%i in (`wsl -d %DISTRO% wslpath "%~dp0provision.sh" 2^>nul`) do set "WSL_SCRIPT=%%i"
if not defined WSL_SCRIPT (
    echo [!] Could not translate script path into WSL. Aborting.
    pause & exit /b 1
)

REM Strip any CR bytes that a Windows editor / git autocrlf may have injected
REM into provision.sh, then pipe the cleaned script into bash. Without this,
REM bash reads `set -euo pipefail\r` and errors on the literal \r.
wsl -d %DISTRO% -u root -- bash -c "sed 's/\r$//' '!WSL_SCRIPT!' | ABACUS_CHINA_MIRROR=!ABACUS_CHINA_MIRROR! bash"
if errorlevel 1 (
    echo [!] Provisioning failed. See output above.
    pause & exit /b 1
)

REM --- 6. Install Windows-side launchers ---
if not exist "%BIN_DIR%" mkdir "%BIN_DIR%"

> "%BIN_DIR%\abacus.cmd" echo @echo off
>> "%BIN_DIR%\abacus.cmd" echo set "WSLENV=OMP_NUM_THREADS:MKL_NUM_THREADS:OPENBLAS_NUM_THREADS:%%WSLENV%%"
>> "%BIN_DIR%\abacus.cmd" echo wsl -d %DISTRO% --cd "%%CD%%" -- abacus %%*

> "%BIN_DIR%\abacus-mpi.cmd" echo @echo off
>> "%BIN_DIR%\abacus-mpi.cmd" echo set "WSLENV=OMP_NUM_THREADS:MKL_NUM_THREADS:OPENBLAS_NUM_THREADS:%%WSLENV%%"
>> "%BIN_DIR%\abacus-mpi.cmd" echo wsl -d %DISTRO% --cd "%%CD%%" -- abacus-mpi %%*

REM Record install state so the uninstaller knows whether we added the distro.
> "%INSTALL_DIR%\install-state.txt" echo distro=%DISTRO%
>> "%INSTALL_DIR%\install-state.txt" echo distro_preexisted=!DISTRO_PREEXISTED!

REM --- 7. Add BIN_DIR to user PATH (idempotent, no 1024-char truncation) ---
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
 "$b='%BIN_DIR%'; $p=[Environment]::GetEnvironmentVariable('PATH','User'); if([string]::IsNullOrEmpty($p)){[Environment]::SetEnvironmentVariable('PATH',$b,'User')} elseif(($p -split ';') -notcontains $b){[Environment]::SetEnvironmentVariable('PATH',$p.TrimEnd(';')+';'+$b,'User')}"

echo.
echo =========================================
echo   Installation complete!
echo =========================================
echo.
echo Open a NEW terminal window, cd into a case directory, then run:
echo     abacus
echo.
echo For parallel execution with N processes:
echo     abacus-mpi -n 4
echo.
echo To enter the Linux shell for advanced use:
echo     wsl -d %DISTRO%
echo.
pause
exit /b 0
