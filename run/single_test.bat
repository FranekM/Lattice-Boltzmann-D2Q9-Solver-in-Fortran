@echo off
setlocal

REM ==========================================================
REM Usage:
REM     run.cmd <num_images>
REM Example:
REM     run.cmd 64
REM ==========================================================

if "%~1"=="" (
    echo Usage: %~nx0 ^<num_images^>
    exit /b 1
)

set N=%~1

set EXE=bin\lbm.exe
set CONFIG=lbm.conf

set FOR_COARRAY_CONFIG_FILE=%CONFIG%

echo.
echo Compiling program...

make -f Makefile.windows clean
make -f Makefile.windows

if errorlevel 1 (
    echo.
    echo Compilation failed.
    exit /b %errorlevel%
)

echo.
echo Running experiment with %N% images

(echo -genvall -genv I_MPI_FABRICS=shm:ofi -n %N% %EXE%) > "%CONFIG%"

%EXE%

if errorlevel 1 (
    echo.
    echo Experiment failed.
    exit /b %errorlevel%
)

echo Experiment completed successfully.
echo.

endlocal