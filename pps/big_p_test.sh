#!/bin/bash
resume_boinc()
# Resume BOINC if it's been stopped.  (Just continues otherwise.)
{
	boinccmd --set_run_mode auto
}

if [ $# = 0 ] ; then echo "Usage: $0 ppsieve_program" ; exit ; fi
if [ "$1" = "all" ] ; then
	sh big_p_test.sh ppsieve-x86_64-linux -s no
	sh big_p_test.sh ppsieve-x86_64-linux -s yes
	sh big_p_test.sh ppsieve-x86-linux-sse2
	sh big_p_test.sh ppsieve-x86-linux-sse2 -s no
	sh big_p_test.sh ppsieve-x86-linux
	exit
fi
if [ ! -x $1 ] ; then echo "Usage: $0 ppsieve_program" ; exit ; fi
if [ -f ppcheck*.txt ] ; then rm ppcheck*.txt ; fi
mv ppfactors.txt ppfactors.txt.bak
if killall boincmgr ; then sleep 1 ; fi
#boinccmd --set_run_mode never
# trap keyboard interrupt (control-c)
#trap resume_boinc SIGINT
sleep 3
./$*  -p80T -P80000050M -ippse_37TE1.txt -z normal -q -t 2 #-fpps_test.txt

resume_boinc
if grep '|' ppfactors.txt | sort | cmp pps_test80t.txt - ; then
	echo Test OK
	mv ppfactors.txt.bak ppfactors.txt
else
	echo TEST INVALID!!!
	printf "See diff (y/n)? "
	read keyin
	if [ $keyin == 'y' ] ; then
		sort ppfactors.txt | diff pps_test80t.txt - | less
	fi
	printf "\nKeep bad ppfactors.txt (y/n)? "
	read keyin
	if [ $keyin != 'y' ] ; then mv ppfactors.txt.bak ppfactors.txt ; fi
fi
