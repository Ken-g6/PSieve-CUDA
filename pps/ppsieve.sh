# Enter your range here:
PMIN=##G
PMAX=##G

# Leave the rest of this stuff alone.
# It automatically finds the best version and number of cores.
SIEVEFILE=ppse_37TE1.txt
PROCS=`grep '^flags' /proc/cpuinfo | wc -l`
MACHINE_TYPE=`uname -m`
if [ "$MACHINE_TYPE" == "x86_64" ]; then
  APP=ppsieve-x86_64-linux
  if grep "AMD" /proc/cpuinfo > /dev/null ; then
    ALT="-a yes"
  fi
else
  if grep '^flags' /proc/cpuinfo | grep sse2 > /dev/null ; then
    APP=ppsieve-x86-linux-sse2
  else
    APP=ppsieve-x86-linux
  fi
fi

if [ "$PMIN" = "##G" -o "$PMAX" = "##G" ] ; then
	echo Please edit this file and enter your range at the top.
	echo Then try again!
	exit
fi

echo Running $APP with $PROCS threads.
if ./$APP -p$PMIN -P$PMAX -i$SIEVEFILE -ffppse_$PMIN-$PMAX.txt -t$PROCS -q $ALT $* ; then
	echo Run completed successfully!
else
	if [ ! -f $SIEVEFILE ] ; then
		echo "Don't forget to download the sieve file!"
		echo "It's at http://pgllr.mine.nu/sieves/ppse/sievefile/"
		exit
	fi

	echo Downloading source
	mkdir src
	cd src
	wget -O src.zip http://home.att.net/~k.brazier/programs/ppsieve-src.zip
	unzip src.zip
	make
	cp pps/$APP ..
	cd ..
	rm -f src/pps/*
	rmdir src/pps
	rm -f src/*
	rmdir src
	if ./$APP -p$PMIN -P$PMAX -ippse_37TE1.txt -ffppse_$PMIN-$PMAX.txt -t$PROCS -q $ALT $* ; then
		echo Run completed successfully!
	else
		echo 'Error after building from source.  Sorry. :('
	fi
fi


