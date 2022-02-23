#!/bin/tcsh

#to be run from $homedir

setenv runid io0030
setenv finyear 1100
setenv ndpyr 360

setenv masterdir $PWD
setenv timmsoutdir $WORK/MisomipPlus/$runid/timmsrun

setenv fesomsourcedir $HOME/MisomipPlus/$runid/fesom
setenv fesommeshdir $WORK/MisomipPlus/$runid/fesommesh
setenv fesomrundir $WORK/MisomipPlus/$runid/fesomrun
setenv fesomdatadir $WORK/MisomipPlus/$runid/fesomdata

setenv uasourcedir $HOME/MisomipPlus/$runid/ua
setenv uadatadir $WORK/MisomipPlus/$runid/uadata
setenv uarstdir $WORK/MisomipPlus/$runid/uarst
setenv uarundir $WORK/MisomipPlus/$runid/uarun

# no users below this line
printf '%s\n' "=====================================================================  ;,//;,    ,;/ ================"
printf '%s\n' "                                  Welcome to Mr. Timms                o:::::::;;///                  "
printf '%s\n' "==================================================================== >::::::::;;\\\ ================="
printf '%s\n' "                                                                       ''\\\\\'' ';\                 "

echo 'masterdir='$masterdir
echo $runid   > $timmsoutdir/$runid.dat

setenv fesomClockFile $fesomdatadir/$runid.clock
setenv timmsClockFile $timmsoutdir/$runid.clock
echo 'FESOM CLOCKFILE is' $fesomClockFile
echo 'Timms CLOCKFILE is' $timmsClockFile

if (-e "$timmsClockFile") then
  echo "======================================="
  echo "extract data from Timms clockfile "
  echo "======================================="
  setenv oldYear `awk -F";" 'NR == 1 { print $1}' $timmsClockFile | awk '{print $1}'`
  setenv oldMonth `awk -F";" 'NR == 1 { print $1}' $timmsClockFile | awk '{print $2}'`
  setenv newYear `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $1}'`
  setenv newMonth `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $2}'`
  setenv yearFrac `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $3}'`

  echo $oldYear
  echo $oldMonth
  echo $newYear
  echo $newMonth
  echo $yearFrac
endif

echo "========================"
echo "go to" $1
echo "========================"
goto $1

start:
ocean2ice:
echo "======================================="
echo "extract data from FESOM clockfile "
echo "======================================="
#VAR11=`awk -F";" 'NR == [ZEILE] { print $1}' $FILE | awk '{print $[SPALTE]}'`
setenv oldTime `awk -F";" 'NR == 1 { print $1}' $fesomClockFile | awk '{print $1}'`
setenv oldDay `awk -F";" 'NR == 1 { print $1}' $fesomClockFile | awk '{print $2}'`
setenv oldYear `awk -F";" 'NR == 1 { print $1}' $fesomClockFile | awk '{print $3}'`
setenv newTime `awk -F";" 'NR == 2 { print $1}' $fesomClockFile | awk '{print $1}'`
setenv newDay `awk -F";" 'NR == 2 { print $1}' $fesomClockFile | awk '{print $2}'`
setenv newYear `awk -F";" 'NR == 2 { print $1}' $fesomClockFile | awk '{print $3}'`
echo "$oldTime"
echo "$oldDay"
echo "$oldYear"
echo "$newTime"
echo "$newDay"
echo "$newYear"

echo "================================================"
echo "update relevant time data in Mr. Timms clockfile"
echo "================================================"

set newMonth = `echo "scale=6; $newDay/$ndpyr*12" | bc | xargs printf "%02.0f"`
if ($newDay == 1) then
    set yearFrac = `printf "%.11f" "0" | cut -d "." -f2`
    set oldMonth = `printf "%02.0f" "11"`
else
    set yearFrac = `echo "scale=12; $newDay / $ndpyr" | bc | xargs printf "%.11f" | cut -d "." -f2 `
    set oldMonth = `echo "scale=6; $newMonth-1" | bc | xargs printf "%02.0f"`
endif
#set oldmonth=`echo "scale=2; $newmonth -1" | bc | xargs printf "%02d"`

echo $oldMonth
echo $newMonth
echo $yearFrac

printf "$oldYear $oldMonth\n$newYear $newMonth $yearFrac" > $timmsClockFile

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after getclockfile has been completed"
 echo "==========================================="
 exit
endif


prepare4ua:
echo "========================"
echo " run prepare4ua"
echo "========================"

echo "preserving old restart file"
cp $uarstdir/RestartMismipPlus-$runid.mat $uarstdir/RestartMismipPlus-$runid.$oldYear.$oldMonth.mat

echo "using "$newYear.$yearFrac " for current ua step."
sed -i "s/CtrlVar.TotalTime=.*/CtrlVar.TotalTime=$newYear.$yearFrac;/" $uasourcedir/DefineInitialInputs.m
sed -i "s/UserVar.StartOutputID=.*/UserVar.StartOutputID='$oldYear.$oldMonth';/" $uasourcedir/DefineInitialInputs.m
sed -i "s/UserVar.FinalOutputID=.*/UserVar.FinalOutputID='$newYear.$newMonth';/" $uasourcedir/DefineInitialInputs.m

setenv fesommeltfile $fesomdatadir/$runid.$oldYear.forcing.diag.nc
setenv fesomcoordfile $fesommeshdir/$oldYear.$oldMonth/nod2d.out
echo "using" $fesommeltfile " and" $fesomcoordfile " for current ua melt rates."
sed -i "s~fesomMeltPath=.*~fesomMeltPath= '$fesommeltfile';~" $uasourcedir/DefineMassBalance.m
sed -i "s~fesomCoordPath=.*~fesomCoordPath= '$fesomcoordfile';~" $uasourcedir/DefineMassBalance.m


if ($status != 0) then
 echo "--------------------------------------------------------------------------------------"
 echo "Attempt to configure DefineInitialInputs.m failed. "
 echo "Exit."
 echo "--------------------------------------------------------------------------------------"
 exit
endif

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after prepare4ua has been completed"
 echo "==========================================="
 exit
endif

launchua:                        
echo "==============="
echo "now launch Ua"
echo "==============="

cd $uasourcedir
sbatch ua.slurm

if ($2 == 'only') then
 echo "======================================================"
 echo "exit after launchua has been completed"
 echo "======================================================"
 exit
endif

setenv newYear `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $1}'`

if ($newYear >= $finyear) then
 echo "==================================================================================="
 echo "exit after launchua has been completed. Final year ("$finyear") has been reached."
 echo "==================================================================================="
 exit
endif


launchualookup:
date
echo "========================================"
echo "launch ua lookup job"
echo "========================================"

setenv oldYear `awk -F";" 'NR == 1 { print $1}' $timmsClockFile | awk '{print $1}'`
setenv oldMonth `awk -F";" 'NR == 1 { print $1}' $timmsClockFile | awk '{print $2}'`
setenv newYear `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $1}'`
setenv newMonth `awk -F";" 'NR == 2 { print $1}' $timmsClockFile | awk '{print $2}'`

cd $masterdir
./check4uadata.csh yes > $timmsoutdir/check4uadata.log &
sleep 2

echo "========================================================================================================="
echo "Ocan2ice part of Mr. Timms script completed; Ua is supposed to run for year $oldYear and month $oldMonth"
echo "========================================================================================================="

exit

# end of ocean to ice
# #######################################################################################################
# # start of ice to ocean    



ice2ocean:
meshgen:
echo "======================"
echo 'call meshgen.m'
echo "======================"

setenv newmeshdir $fesommeshdir/$newYear.$newMonth
setenv uaresultfile $uadatadir/$newYear.$newMonth-Nodes*.mat
setenv newgoodfile $fesommeshdir/meshgen.goodfile.$newYear.$newMonth
setenv oldgoodfile $fesommeshdir/meshgen.goodfile.$oldYear.$oldMonth

sed -i "s~meshOutPath=.*~meshOutPath='$newmeshdir/';~" meshgen.m
sed -i "s~Ua_path=.*~Ua_path='$uaresultfile';~" meshgen.m
sed -i "s~goodfile_path=.*~goodfile_path='$newgoodfile';~" meshgen.m
mkdir -p $newmeshdir/dist

matlab.sh -s -S"-wprod-0304" -M"-nojvm -r run('meshgen.m')"  #> meshgen.$yearfromicemod.log
 #matlab.sh -s -M" -nodisplay -r run('~/test.m')" -S"--time=6:00:00 -c4 --mem=10000"

if (-e $newgoodfile) then
 echo "||||||||||||||||||||||||||||"
 echo "returned from meshgen.m ok; remove old goodfile"
 echo "||||||||||||||||||||||||||||"
 rm $oldgoodfile
else
 echo ----------------------------------------------
 echo 'returned from meshgen.m accidentally: stop'
 echo ----------------------------------------------
 exit
endif

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after meshgen has been completed"
 echo "==========================================="
 exit
endif

remapfesom:
echo "================================================="
echo 'now remap fesom data to new mesh'
echo "================================================="

setenv newmeshdir $fesommeshdir/$newYear.$newMonth
setenv oldmeshdir $fesommeshdir/$oldYear.$oldMonth

setenv oldocefile $fesomdatadir/$runid.$oldYear.oce.nc
setenv oldicefile $fesomdatadir/$runid.$oldYear.ice.nc
setenv newocefile $fesomdatadir/$runid.$oldYear.oce.ini.nc
setenv newicefile $fesomdatadir/$runid.$oldYear.ice.ini.nc

echo "from: "
echo $oldmeshdir
echo $oldocefile
echo $oldicefile

echo ""
echo "to: "
echo $newmeshdir
echo $newocefile
echo $newicefile

sed -i "s~oldOceFile=.*~oldOceFile='$oldocefile';~" remap.m
sed -i "s~newOceFile=.*~newOceFile='$newocefile';~" remap.m
sed -i "s~oldIceFile=.*~oldIceFile='$oldicefile';~" remap.m
sed -i "s~newIceFile=.*~newIceFile='$newicefile';~" remap.m
sed -i "s~oldMeshPath=.*~oldMeshPath='$oldmeshdir/';~" remap.m
sed -i "s~newMeshPath=.*~newMeshPath='$newmeshdir/';~" remap.m

matlab.sh -s -S"-wprod-0304" -M"-nojvm -r run('remap.m')"

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after remap_data has been completed"
 echo "==========================================="
 exit
endif

prepare4fesom:
echo "============="
echo 'prepare fesom'
echo "============="

echo "archive fesom output using year.month format ... "
set names = ('forcing.diag' 'ice.diag' 'ice.mean' 'mesh.diag' 'oce.diag' 'oce.mean' 'oce' 'ice')
foreach name ($names)
   mv $fesomdatadir/$runid.$oldYear.$name.nc $fesomdatadir/$runid.$oldYear.$oldMonth.$name.nc
   echo "$name is done"
end
set names = ('ice.ini' 'oce.ini')
foreach name ($names)
   cp $fesomdatadir/$runid.$oldYear.$name.nc $fesomdatadir/$runid.$oldYear.$oldMonth.$name.nc
   echo "$name is done"
end
sleep 5

setenv newmeshdir $fesommeshdir/$newYear.$newMonth
setenv configfile $fesomsourcedir/namelist.config
echo 'now modify namelist.config using year '$oldYear' and mesh '$newmeshdir
#cp $configfile $fesomrundir/namelist.config.$yeartoicemod
sed -i "s~yearnew=.*~yearnew=$oldYear~" $configfile
sed -i "s~MeshPath=.*~MeshPath='$newmeshdir/'~" $configfile

echo 'remove mesh/dist folder from last year'
rm -r $fesommeshdir/$oldYear.$oldMonth/dist

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after prepare_fesom has been completed"
 echo "==========================================="
 exit
endif

launchfesom:
echo "=========================================="
echo 'now launch fesom'
echo "=========================================="

cd $fesomsourcedir
sbatch fesom.slurm | cut -d ' ' -f4 | tee $timmsoutdir/$runid.jobid.dat

if ($status != 0) then
 echo "-----------------------------------------"
 echo "Attempt to launch FESOM has failed. Exit."
 echo "-----------------------------------------"
 exit
endif

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after launch_fesom has been completed"
 echo "==========================================="
 exit
endif

launchfesomlookup:
date
echo "=========================================="
echo 'now launch sleeping lookup job'
echo "=========================================="

cd $masterdir
./check4fesomdata.csh yes > $timmsoutdir/check4fesomdata.log &
sleep 2

echo "================================================="
echo "Mr. Timms script completed; FESOM now runs for year $newYear and month $newMonth"
echo "================================================="

exit

