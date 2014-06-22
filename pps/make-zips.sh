appname=ppsieve-cl
apppath=pps
rm -f $appname.zip $appname-src.zip
# Make sure everything's committed.
status=`git status | grep nothing | wc -l`
if [ $status == "0" ] ; then
	git status
	exit
fi
# Save CHANGES.txt
git log | grep "^\(    \|Date:\)" > CHANGES.txt
zip -9 $appname.zip $appname-x86-linux $appname-boinc-x86-linux $appname-x86_64-linux $appname-boinc-x86_64-linux $appname-boinc-x86-windows.exe $appname-x86-windows.exe $appname-boinc-x86_64-windows.exe $appname-x86_64-windows.exe README.txt CHANGES.txt $appname.sh $appname.bat *config.txt license.txt
#cd ..
#zip -9 $apppath/$appname-src.zip Makefile $apppath/Makefile *.txt *.[ch] *.cpp *.sln $apppath/*.[ch] $apppath/*.cu $apppath/README.txt $apppath/CHANGES.txt $apppath/*config.txt $apppath/make-*.sh $apppath/$appname.sh $apppath/*.bat vc/*
