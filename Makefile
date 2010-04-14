all: do_pps do_pps_boinc

do_pps:
	cd pps; bash ./make-bins.sh

do_pps_boinc:
	cd pps; bash ./make-bins.sh boinc

clean:
	rm pps/ppsieve-boinc-x86*-linux pps/ppsieve-x86*-linux pps/ppsieve-x86*-windows.exe
