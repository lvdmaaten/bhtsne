git clone https://github.com/hpicgs/bhtsne bhtsne-performance-test --branch dev
cd bhtsne-performance-test
mkdir build
cd build
cmake ..
cmake --build .
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Enterprise\Common7\Tools\VsDevCmd.bat"
cd %~dp0\bhtsne-performance-test\build
devenv /Build Release bhtsne.sln
xcopy ..\..\data\* Release\
cd Release
performance_test.exe
xcopy performance_*.csv ..\..\..\results\
pause