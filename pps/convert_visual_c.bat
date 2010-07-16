move ..\main.c ppsieve-cuda.c
copy ..\*.c .
copy ..\*.h .
copy ..\vc\*.* .
cd ..
start ppsieve-cuda.sln
