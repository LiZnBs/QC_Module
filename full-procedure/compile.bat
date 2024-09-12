@echo off
wolframscript -file .\Input.wls
cd .\Real-DSE
REM compile wolfram files: Real-DSE
wolframscript -file .\Main.wls
move /Y .\Result\Real-dse_A-B_p2_All-results.wls ..\Complex-DSE\Parameters\Real-dse_A-B_p2_All-results.wls
move /Y .\Result\Real-dse_Z2-Z4.wls ..\Complex-DSE\Parameters\Real-dse_Z2-Z4.wls
move /Y .\Result\Z2.txt ..\BSE\DSE_Results\Z2.txt
cd ..\Complex-DSE
REM compile wolfram files: Complex-DSE
wolframscript -file .\Main.wls
move /Y .\Result\Complex-dse_A.txt ..\BSE\DSE_Results\Complex-dse_A.txt
move /Y .\Result\Complex-dse_B.txt ..\BSE\DSE_Results\Complex-dse_B.txt
move /Y .\Result\kp.txt ..\BSE\DSE_Results\kp.txt
move /Y .\Result\zp.txt ..\BSE\DSE_Results\zp.txt
move /Y .\Result\yp.txt ..\BSE\DSE_Results\yp.txt
cd ..\BSE
REM compile fortran files: BSE
gfortran -fopenmp Main.f95 -o Main
.\Main