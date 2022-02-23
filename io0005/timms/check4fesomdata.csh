#!/bin/csh

setenv fesomdatadir $1
setenv runid $2
setenv yearfromicemod $3
setenv firsttime $4
setenv timmsoutdir $5

date

if ($firsttime == yes) then
  echo "firsttime=yes: sleep for 10 seconds"
  sleep 10 
  #sleep 2
else
  echo "firsttime=no: sleep for 10 seconds"
  sleep 10
  #sleep 1 
endif

set ID=`cat $timmsoutdir/$runid.jobid.dat`
set jobfile="$fesomrundir/slurm-$ID.out"

if (-e "$jobfile") then
  echo "Checking for 'not converged' or 'Velocity too large' statement in $jobfile."
  if ( `grep -c "not converged" $jobfile` != 0 || `grep -c "Velocity too large" $jobfile` != 0) then
    echo "DETECTED! CANCEL CURRENT SLURM JOB."
    scancel $ID
    mv $jobfile $jobfile.canceled
    set value=`grep "step_per_day" $fesomrundir/namelist.config | cut -d '=' -f2 | cut -d '!' -f1`
    if ( $value < 960 ) then
       @ newValue= $value * 2
       echo "adjusting step_per_day from $value to $newValue."
       sed -i "s~step_per_day=.*~step_per_day=$newValue~" $fesomrundir/namelist.config
       echo "Relaunch fesom with adjusted time step."
       cd $fesomrundir
       sbatch fesom.slurm | cut -d ' ' -f4 | tee $timmsoutdir/$runid.jobid.dat
       echo "check again in 1 minutes"
       cd $masterdir
       ./check4fesomdata.csh $fesomdatadir $runid $yearfromicemod no $timmsoutdir >> $timmsoutdir/check4fesomdata.log &
       exit
    else
       echo "Time step_per_day is already >= 960. We do nothing and stop the loop."
       exit
    endif
  else
    echo "No 'not coverged' statement found. Continue running job $ID and check if run has finished successfully."
  endif
else
   echo "$jobfile does not yet exist. Check again in 1 minute."
   ./check4fesomdata.csh $fesomdatadir $runid $yearfromicemod no $timmsoutdir >> $timmsoutdir/check4fesomdata.log &
   exit
endif

#Now check if run has finished.
setenv CLOCKFILE $runid.clock
@ day2wait4 = (`echo $yearfromicemod | cut -d "." -f2` + 10) * 360 / 100
@ year2wait4 = `echo $yearfromicemod | cut -d "." -f1` + $day2wait4 / 360
echo "The day to wait for is " $day2wait4
echo "The year to wait for is " $year2wait4

echo "now get " $fesomdatadir/$CLOCKFILE
setenv OLDDAY `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
setenv OLDYEAR `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
setenv NEWDAY `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
setenv NEWYEAR `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
echo "OLDDAY and NEWYEAR now are" $OLDDAY $NEWYEAR
if ( $NEWYEAR == $year2wait4 && $OLDDAY == $day2wait4 ) then
   set value=`grep "step_per_day" $fesomrundir/namelist.config | cut -d '=' -f2 | cut -d '!' -f1`
   if ( $value > 240 ) then
       @ newValue= $value / 2 
       echo "adjusting step_per_day from $value to $newValue."
       sed -i "s~step_per_day=.*~step_per_day=$newValue~" $fesomrundir/namelist.config
   else
       echo "Time step_per_day is already = 120."
   endif
   echo "launch timms.csh start all"
   ./timms.csh start all > $timmsoutdir/timms_oce2ice_$yearfromicemod.log &
   echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
   echo "Mr. Timms oce2ice has been launched: exit check4fesomdata loop"
   echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
   exit
else
   echo "Check again in 1 minute."
   ./check4fesomdata.csh $fesomdatadir $runid $yearfromicemod no $timmsoutdir >> $timmsoutdir/check4fesomdata.log &
   echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n "
endif
exit

