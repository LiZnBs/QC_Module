@echo off
goto main
REM functions
REM excecuting normailzation
:normroutine
    REM for excecuting normalization
    cd .\normalization
    copy /Y ..\BSE\DSE_Results\Complex-dse_A.txt .\DSE-BSE_results\Complex-dse_A1.txt
    copy /Y ..\BSE\DSE_Results\Complex-dse_B.txt .\DSE-BSE_results\Complex-dse_B1.txt
    move /Y ..\BSE\Complex-f.txt .\DSE-BSE_results\Complex-f.txt
    REM excecuting adjacent M complex-DSE
    cd ..\Complex-DSE
    wolframscript -file .\decay-constant-api.wls
    wolframscript -file .\Main.wls
    move /Y .\Result\Complex-dse_A.txt ..\normalization\DSE-BSE_results\Complex-dse_A2.txt
    move /Y .\Result\Complex-dse_B.txt ..\normalization\DSE-BSE_results\Complex-dse_B2.txt
    REM continue excecuting normalization
    cd ..\normalization
    .\bin\main.exe
    copy /Y .\results\normalized-f.txt ..\lapton-decay-constant\DSE-BSE_results\normalized-f.txt
    copy /Y .\results\normalized-f.txt ..\rp_vacuum-quark-condensate\DSE-BSE_results\normalized-f.txt
    set flag=1
goto :eof

REM main
:main
cd %~dp0
wolframscript -file .\Input.wls
cd .\Real-DSE
REM compile wolfram files: Real-DSE
wolframscript -file .\Main.wls
move /Y .\Result\Real-dse_A-B_p2_All-results.wls ..\Complex-DSE\Parameters\Real-dse_A-B_p2_All-results.wls
move /Y .\Result\Real-dse_Z2-Z4.wls ..\Complex-DSE\Parameters\Real-dse_Z2-Z4.wls
copy /Y .\Result\Z2.txt ..\BSE\DSE_Results\Z2.txt
cd ..\Complex-DSE
REM compile wolfram files: Complex-DSE
wolframscript -file .\Main.wls
move /Y .\Result\Complex-dse_A.txt ..\BSE\DSE_Results\Complex-dse_A.txt
move /Y .\Result\Complex-dse_B.txt ..\BSE\DSE_Results\Complex-dse_B.txt
cd ..\BSE
REM compile fortran files: BSE
gfortran -fopenmp Main.f95 -o Main
.\Main
REM flag use to check whether have calculated normalization or not 
cd %~dp0
set flag=0 
:input1
set /p userinput1=Whether do you want to calculate the lapton decay constant? (y/n):
if "%userinput1%"=="y" (
    call :normroutine
    REM continue excecuting decay constant
    cd ..\lapton-decay-constant
    copy /Y ..\BSE\DSE_Results\Z2.txt .\DSE-BSE_results\Z2.txt
    copy /Y ..\BSE\DSE_Results\Complex-dse_A.txt .\DSE-BSE_results\Complex-dse_A.txt
    copy /Y ..\BSE\DSE_Results\Complex-dse_B.txt .\DSE-BSE_results\Complex-dse_B.txt
    .\bin\main.exe
) else (
    if "%userinput1%"=="n" (
        REM exit
        
        ) else (
    echo Invalid input
    goto input1
    )
)
cd %~dp0
:input2
set /p userinput2=Whether do you want to calculate the rp (vacuum-quark-condensate)? (y/n):
if "%userinput2%"=="y" (
    if %flag%==0 (
        call :normroutine
    )
    REM excecuting rp
    cd %~dp0
    cd .\rp_vacuum-quark-condensate
    copy /Y ..\Real-DSE\Result\Z4.txt .\DSE-BSE_results\Z4.txt
    copy /Y ..\BSE\DSE_Results\Complex-dse_A.txt .\DSE-BSE_results\Complex-dse_A.txt
    copy /Y ..\BSE\DSE_Results\Complex-dse_B.txt .\DSE-BSE_results\Complex-dse_B.txt
    call .\compile_rp.bat
) else (
    if "%userinput2%"=="n" (
        REM exit

        ) else (
    echo Invalid input
    goto input2
    )
)
