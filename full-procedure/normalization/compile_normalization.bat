@echo off
cd %~dp0
REM set the directories
set SOURCE_DIR=.\source
set BIN_DIR=.\bin
set MOD_DIR=.\mod
set OBJ_DIR=.\obj

REM Create bin directory if not exist
if not exist bin mkdir bin 
if not exist mod mkdir mod
if not exist obj mkdir obj

REM compile fortran files
gfortran -fopenmp -J %MOD_DIR% -c %SOURCE_DIR%\parameters.f95 -o %OBJ_DIR%\parameters.o
gfortran -fopenmp -J %MOD_DIR% -c %SOURCE_DIR%\functions.f95 -o %OBJ_DIR%\functions.o
gfortran -fopenmp -J %MOD_DIR% -c %SOURCE_DIR%\main.f95 -o %OBJ_DIR%\main.o
REM link the object files
gfortran -fopenmp -o %BIN_DIR%\main %OBJ_DIR%\parameters.o %OBJ_DIR%\functions.o %OBJ_DIR%\main.o
REM execute the program if the compilation was successful
if %ERRORLEVEL% EQU 0 .\bin\main