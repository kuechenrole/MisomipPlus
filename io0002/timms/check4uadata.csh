#!/bin/csh

setenv uadatadir $1
setenv year2wait4 $2
setenv firsttime $3
setenv timmsoutdir $4

date

if ($firsttime == yes) then
  echo "firsttime=yes: sleep for 30 seconds"
  sleep 30
else
  echo "firsttime=no: sleep for 10 seconds"
  sleep 10 
endif

date
echo "The year to wait for is " $year2wait4

ls ${uadatadir}/0${year2wait4}00-Nodes*.mat
if ($status == 0) then
   sleep 2
   echo "Ua results file exists: launch Mr. Timms ice2ocean"
   ./timms.csh ice2ocean all > $timmsoutdir/timms_ice2oce_$year2wait4.log &

   echo "setting UserVar.CouplingStart=0."
   sed -i "s~UserVar.CouplingStart=.*~UserVar.CouplingStart=0;~" $uasourcedir/DefineInitialInputs.m

   echo "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
   echo "MrTimms.csh ice2ocean has been launched, check4uadata terminates succefully"
   echo "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
else
   echo "Ua results does not exist: check again in 1 minute."
   ./check4uadata.csh $uadatadir $year2wait4 no $timmsoutdir >> $timmsoutdir/check4uadata.log &
endif
exit

