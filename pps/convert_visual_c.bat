move ..\main.c ppsieve-cuda.c
copy ..\*.c .
copy ..\*.cpp .
copy ..\*.h .
copy ..\vc\*.* .
cd ..
start ppsieve-cuda.sln
