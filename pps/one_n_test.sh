#!/bin/bash
resume_boinc()
# Resume BOINC if it's been stopped.  (Just continues otherwise.)
{
	boinccmd --set_run_mode auto
}

if [ $# = 0 ] ; then echo "Usage: $0 ppsieve_program [args]" ; exit ; fi
if [ "$1" = "all" ] ; then
	sh one_n_test.sh ppsieve-x86_64-linux $2 $3 -s yes
	sh one_n_test.sh ppsieve-x86_64-linux $2 $3 -s no
	sh one_n_test.sh ppsieve-x86-linux-sse2 $2 $3
	sh one_n_test.sh ppsieve-x86-linux-sse2 $2 $3 -s no
	sh one_n_test.sh ppsieve-x86-linux $2 $3
	exit
fi
if [ ! -x $1 ] ; then echo "Usage: $0 ppsieve_program" ; exit ; fi
if [ -f ppcheck*.txt ] ; then rm ppcheck*.txt ; fi
mv ppfactors.txt ppfactors.txt.bak
if killall boincmgr ; then sleep 1 ; fi
boinccmd --set_run_mode never
# trap keyboard interrupt (control-c)
#trap resume_boinc SIGINT
sleep 3
tset=''
if [ "$2" = "-t" -o "$3" = "-t" -o "$4" = "-t" -o "$5" = "-t" ] ; then
	tset=''
else
	tset='taskset 0x00000002'
fi

$tset ./$*  -p447000000 -P447000100 -k1200 -K9999 -n2M -N3M -z normal -q

resume_boinc
if grep '|' ppfactors.txt | sort | cmp oneOKsort.txt - ; then
	echo Test OK
	mv ppfactors.txt.bak ppfactors.txt
else
	echo TEST INVALID!!!
	printf "See diff (y/n)? "
	read keyin
	if [ $keyin == 'y' ] ; then
		sort ppfactors.txt | gvimdiff -- oneOKsort.txt -
	fi
	printf "\nKeep bad ppfactors.txt (y/n)? "
	read keyin
	if [ $keyin != 'y' ] ; then mv ppfactors.txt.bak ppfactors.txt ; fi
fi
