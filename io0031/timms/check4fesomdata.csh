#!/bin/csh

setenv firsttime $3

date

if ($firsttime == yes) then
  echo "firsttime=yes: sleep for 20 seconds"
  sleep 20 
  #sleep 2
else
  echo "firsttime=no: sleep for 10 seconds"
  sleep 10
  #sleep 1 
endif

set ID=`cat $timmsoutdir/$runid.jobid.dat`
set jobfile="$fesomrundir/slurm-$ID.out"
set configfile = "$fesomsourcedir/namelist.config"

if (! -e "$jobfile") then
   echo "$jobfile does not yet exist. Check again in 1 minute."
   ./check4fesomdata.csh no >> $timmsoutdir/check4fesomdata.log &
   exit
else
  echo "Checking for 'not converged' or 'Velocity too large' statement in $jobfile."
  if ( `grep -c "not converged" $jobfile` != 0 || `grep -c "Velocity too large" $jobfile` != 0) then
    echo "DETECTED! CANCEL CURRENT SLURM JOB."
    scancel $ID
    mv $jobfile $jobfile.canceled
    set value=`grep "step_per_day" $configfile | cut -d '=' -f2 | cut -d '!' -f1`
    if ( $value < 2880 ) then
       @ newValue= $value * 2
       echo "adjusting step_per_day from $value to $newValue."
       sed -i "s~step_per_day=.*~step_per_day=$newValue~" $configfile
       echo "Relaunch fesom with adjusted time step."
       cd $fesomsourcedir
       sbatch fesom.slurm | cut -d ' ' -f4 | tee $timmsoutdir/$runid.jobid.dat
       echo "check again in 1 minutes"
       cd $masterdir
       ./check4fesomdata.csh no >> $timmsoutdir/check4fesomdata.log &
       exit
    else
       echo "Time step_per_day is already >= 2880. We do nothing and stop the loop."
       exit
    endif
  else
    echo "No 'not coverged' statement found. Continue running job $ID and check if run has finished successfully."
    if ( `grep -c "successfully completed" $jobfile` != 0 ) then
      echo "FESOM has finished."
      set value=`grep "step_per_day" $configfile | cut -d '=' -f2 | cut -d '!' -f1`
      if ( $value > 240 ) then
        @ newValue= $value / 2 
        echo "adjusting step_per_day from $value to $newValue."
        sed -i "s~step_per_day=.*~step_per_day=$newValue~" $configfile
      else
        echo "Time step_per_day is already = 240."
      endif
      echo "launch timms.csh start all"
      ./timms.csh start all > $timmsoutdir/timms_oce2ice_$newYear.$newMonth.log &
      echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
      echo "Mr. Timms oce2ice has been launched: exit check4fesomdata loop"
      echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
      exit
    else
      echo "No 'successfully completed' found yet. Check again in 1 minute."
      ./check4fesomdata.csh no >> $timmsoutdir/check4fesomdata.log &
      echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n "
    endif
  endif
endif
exit


