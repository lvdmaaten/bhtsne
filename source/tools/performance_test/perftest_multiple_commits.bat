@echo off
if not exist commits.txt (
	echo Please provide a file named commits.txt which contains the commits which should be tested.
	pause
	exit 1
)
if not exist data (
	echo Please provide a data folder containing the test datasets.
	pause
	exit 1
)
echo Started test - %time%
echo|set /p="Cloning repository... "
git clone https://github.com/hpicgs/bhtsne bhtsne-performance-test --branch dev -q
echo Done
echo|set /p="Setting up VS env... "
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Enterprise\Common7\Tools\VsDevCmd.bat" > nul
echo Done
cd %~dp0\bhtsne-performance-test\
for /F "tokens=*" %%A in (..\commits.txt) do (
	echo Processing commit %%A:
	git checkout %%A -q
	mkdir build
	cd build
	echo|set /p="Running cmake... "
	cmake .. > nul
	cmake --build . > nul
	echo Done
	echo|set /p="Compiling... "
	devenv /Build Release bhtsne.sln > nul
	echo Done
	xcopy ..\..\data\* Release\data\ > nul
	cd Release
	echo|set /p="Running tests... "
	performance_test.exe 0 1 %%A
	echo Done
	xcopy performance_*.csv ..\..\..\results\ > nul
	cd ..\..
	rmdir build /s /q
	git clean -d -f -q
	echo Finished processing commit %%A!
)
cd ..
echo Cleaning up...
rmdir bhtsne-performance-test /s /q
echo Done - %time%
pause
