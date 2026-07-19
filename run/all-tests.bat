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
set CONFIG=run\lbm.conf
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


for %%P in (
    shear-wave-decay
	density-wave-decay
	gaussian-wave-decay
	couette-flow
	poiseuille-flow
	sliding-lid
	moving-walls
) do (

    echo.
    echo Running experiment: %%P
    echo Images: %N%

    (echo -genvall -genv I_MPI_FABRICS=shm:ofi -n %N% %EXE% %%P) > "%CONFIG%"

    %EXE% %%P

    if errorlevel 1 (
        echo.
        echo Experiment %%P failed.
    )
)

echo All experiments completed successfully.
echo.

endlocal