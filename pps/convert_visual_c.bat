move ..\main.c ppsieve-cuda.c
copy ..\*.c .
copy ..\*.cpp .
copy ..\*.h .
copy ..\vc\*.* .
perl cltoh.pl appcl.cl > temp.h
if errorlevel 1 goto perlerr
copy appcl_start.h+temp.h > appcl.h
cd ..
start ppsieve-cuda.sln
start ppsieve-cl.sln
goto end

:perlerr
echo Error!  Please install a Perl interpreter to continue.
echo Try http://strawberryperl.com/
:end
