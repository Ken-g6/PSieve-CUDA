all: do_pps

do_pps:
	cd pps; sh ./make-bins.sh

clean:
	rm pps/ppsieve-x86*-linux pps/ppsieve-x86*-windows.exe
