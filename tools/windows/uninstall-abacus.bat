@echo off
setlocal EnableExtensions EnableDelayedExpansion
title ABACUS Uninstaller

set "DISTRO=Ubuntu-22.04"
set "INSTALL_DIR=%LOCALAPPDATA%\ABACUS"
set "BIN_DIR=%INSTALL_DIR%\bin"
set "STATE_FILE=%INSTALL_DIR%\install-state.txt"

REM --- Read state file to learn whether the installer added the distro ---
set "DISTRO_PREEXISTED=1"
if exist "%STATE_FILE%" (
    for /f "usebackq tokens=1,2 delims==" %%a in ("%STATE_FILE%") do (
        if /i "%%a"=="distro_preexisted" set "DISTRO_PREEXISTED=%%b"
    )
)

echo This will remove:
echo     - Launchers under %BIN_DIR%
echo     - %BIN_DIR% from your user PATH
if "!DISTRO_PREEXISTED!"=="0" (
    echo.
    echo The installer added the '%DISTRO%' WSL distribution.
    echo You can either remove the whole distribution, or just remove the
    echo ABACUS files inside it and keep the Linux environment.
) else (
    echo     - /opt/abacus-miniforge inside WSL ^(%DISTRO%^)
    echo.
    echo The %DISTRO% distribution existed before installation and will NOT be removed.
)
echo.

set /p CONFIRM="Continue? (y/N): "
if /i not "!CONFIRM!"=="y" (
    echo Aborted.
    exit /b 0
)

set "REMOVED_DISTRO=0"
if "!DISTRO_PREEXISTED!"=="0" (
    echo.
    set /p REMOVE_DISTRO="Remove the entire '%DISTRO%' WSL distribution? This wipes all files inside it. (y/N): "
    if /i "!REMOVE_DISTRO!"=="y" (
        echo [*] Unregistering %DISTRO%...
        wsl --unregister %DISTRO%
        if not errorlevel 1 set "REMOVED_DISTRO=1"
    )
)

if "!REMOVED_DISTRO!"=="0" (
    where wsl >nul 2>&1
    if not errorlevel 1 (
        wsl -d %DISTRO% -u root -- /bin/true >nul 2>&1
        if not errorlevel 1 (
            echo [*] Removing conda env and launchers inside %DISTRO%...
            wsl -d %DISTRO% -u root -- bash -c "rm -rf /opt/abacus-miniforge /usr/local/bin/abacus /usr/local/bin/abacus-mpi"
        )
    )
)

if exist "%BIN_DIR%"     rd /s /q "%BIN_DIR%"
if exist "%INSTALL_DIR%" rd /s /q "%INSTALL_DIR%"

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
 "$b='%BIN_DIR%'; $p=[Environment]::GetEnvironmentVariable('PATH','User'); if($p){$n=(($p -split ';') | Where-Object { $_ -and ($_ -ne $b) }) -join ';'; [Environment]::SetEnvironmentVariable('PATH',$n,'User')}"

echo.
echo Uninstall complete.
pause
exit /b 0
