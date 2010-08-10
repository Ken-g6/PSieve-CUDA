all:
	cd pps; make

do_pps:
	cd pps; make non-boinc

do_pps_boinc:
	cd pps; make boinc

clean:
	cd pps; make clean
