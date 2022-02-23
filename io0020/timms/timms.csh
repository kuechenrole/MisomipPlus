#!/bin/tcsh

#to be run from $homedir

setenv runid io0020
setenv finyear 1100

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
#cd $masterdir
echo $runid   > $timmsoutdir/$runid.dat

setenv CLOCKFILE $runid.clock
echo 'defined CLOCKFILE to be' $fesomdatadir/$CLOCKFILE


echo "========================"
echo "go to" $1
echo "========================"
goto $1


#Seiteneinstiege:                            
#getclockfile:
#echo "========================"
#echo "getclockfile"
#echo "========================"

#echo 'get file from' $fesomdatadir
#cp $fesomdatadir/$CLOCKFILE .
#if ($status != 0) then
# echo "--------------------------------------------------------------------------------------"
# echo "Attempt to copy" $fesomdatadir/$CLOCKFILE " to masterdir failed."
# echo "Exit."
# echo "--------------------------------------------------------------------------------------"
# exit
#endif

start:
ocean2ice:
echo "======================================="
echo "extract data from clockfile "
echo "======================================="
#VAR11=`awk -F";" 'NR == [ZEILE] { print $1}' $FILE | awk '{print $[SPALTE]}'`
setenv OLDTIME `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $1}'`
setenv OLDDAY `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
setenv OLDYEAR `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
setenv NEWTIME `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $1}'`
setenv NEWDAY `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
setenv NEWYEAR `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
echo "$OLDTIME"
echo "$OLDDAY"
echo "$OLDYEAR"
echo "$NEWTIME"
echo "$NEWDAY"
echo "$NEWYEAR"
@ yearfrac = $OLDDAY * 100 / 360
setenv NEWYEAR $NEWYEAR.`printf "%02d" $yearfrac | tail -c 2`
@ yearfrac = $yearfrac - 10
setenv OLDYEAR $OLDYEAR.`printf "%02d" $yearfrac | tail -c 2`
echo "$OLDYEAR"
echo "$NEWYEAR"
#if ($NEWDAY != 1 || $OLDYEAR == $NEWYEAR || ($OLDDAY != 360 && $OLDDAY != 365)) then
# echo "--------------------------------------------"
# echo "ERROR: invalid clock file - Mr. Timms terminates"
# echo "--------------------------------------------"
#endif

echo $OLDYEAR > $timmsoutdir/$runid.yeartoicemod.dat
echo $NEWYEAR > $timmsoutdir/$runid.yearfromicemod.dat

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


setenv OLDYEAR `cat $timmsoutdir/$runid.yeartoicemod.dat`
setenv NEWYEAR `cat $timmsoutdir/$runid.yearfromicemod.dat`

#echo "Preserving log file"
#cp $uadatadir/ice.log $uadatadir/ice.$OLDYEAR.log
echo "preserving old restart file"
#set uaid=`grep "UserVar.MisExperiment=" $uadatadir/DefineInitialInputs.m | cut -d "'" -f2`
cp $uarstdir/RestartMismipPlus-$runid.mat $uarstdir/RestartMismipPlus-$runid.$OLDYEAR.mat

echo "using years" $OLDYEAR" and" $NEWYEAR" for current ua step."
sed -i "s/CtrlVar.TotalTime=.*/CtrlVar.TotalTime=$NEWYEAR;/" $uasourcedir/DefineInitialInputs.m
#sed -i "s/CtrlVar.RestartTime=.*/CtrlVar.RestartTime=$OLDYEAR;/" $uadatadir/DefineInitialInputs.m


setenv fesommeltfile $fesomdatadir/$runid.`echo $OLDYEAR | cut -d "." -f1`.forcing.diag.nc
setenv fesomcoordfile $fesommeshdir/$OLDYEAR/nod2d.out
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

launchua:                          # !!does not work as an entry on AWI!!
echo "==============="
echo "now launch Ua"
echo "==============="

setenv OLDYEAR `cat $timmsoutdir/$runid.yeartoicemod.dat`
setenv NEWYEAR `cat $timmsoutdir/$runid.yearfromicemod.dat`

cd $uasourcedir
sbatch ua.slurm

if ($2 == 'only') then
 echo "======================================================"
 echo "exit after launchua has been completed"
 echo "======================================================"
 exit
endif

if (`echo $NEWYEAR | cut -d "." -f1` >= $finyear) then
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
cd $masterdir
setenv yearfromicemod `cat $timmsoutdir/$runid.yearfromicemod.dat`

./check4uadata.csh $uadatadir $yearfromicemod yes $timmsoutdir > $timmsoutdir/check4uadata.log &
sleep 2

echo "=========================================================================="
echo "Ocan2ice part of Mr. Timms script completed; Ua is supposed to run for year"
cat $timmsoutdir/$runid.yearfromicemod.dat
echo "=========================================================================="

exit

# end of ocean to ice
# #######################################################################################################
# # start of ice to ocean    



ice2ocean:
meshgen:
echo "======================"
echo 'call meshgen.m'
echo "======================"
setenv yearfromicemod `cat $timmsoutdir/$runid.yearfromicemod.dat`
setenv newmeshdir $fesommeshdir/$yearfromicemod
setenv uaresultfile $uadatadir/$yearfromicemod-Nodes*.mat
setenv goodfile $fesommeshdir/meshgen.goodfile.$yearfromicemod

sed -i "s~meshOutPath=.*~meshOutPath='$newmeshdir/';~" meshgen.m
sed -i "s~Ua_path=.*~Ua_path='$uaresultfile';~" meshgen.m
sed -i "s~goodfile_path=.*~goodfile_path='$goodfile';~" meshgen.m
mkdir -p $newmeshdir/dist
matlab.sh -s -S"-wprod-0304" -M"-nojvm -r run('meshgen.m')"  #> meshgen.$yearfromicemod.log
 #matlab.sh -s -M" -nodisplay -r run('~/test.m')" -S"--time=6:00:00 -c4 --mem=10000"

if (-e $goodfile) then
 echo "||||||||||||||||||||||||||||"
 echo "returned from meshgen.m ok"
 echo "||||||||||||||||||||||||||||"
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

remap_data:
echo =================================================
echo 'now remap fesom data to new mesh'
echo =================================================
setenv yeartoicemod `cat $timmsoutdir/$runid.yeartoicemod.dat`
setenv yearfromicemod `cat $timmsoutdir/$runid.yearfromicemod.dat`
echo 'yeartoicemod=' $yeartoicemod
echo 'yearfromicemod=' $yearfromicemod

setenv newmeshdir $fesommeshdir/$yearfromicemod
setenv oldmeshdir $fesommeshdir/$yeartoicemod

#mkdir $fesomdatadir/arch
#setenv archocefile $fesomdatadir/arch/$runid.$yeartoicemod.oce.nc
#setenv archicefile $fesomdatadir/arch/$runid.$yeartoicemod.ice.nc
setenv oldocefile $fesomdatadir/$runid.`echo $yeartoicemod | cut -d "." -f1`.oce.nc
setenv oldicefile $fesomdatadir/$runid.`echo $yeartoicemod | cut -d "." -f1`.ice.nc
setenv newocefile $fesomdatadir/$runid.`echo $yeartoicemod | cut -d "." -f1`.oce.ini.nc
setenv newicefile $fesomdatadir/$runid.`echo $yeartoicemod | cut -d "." -f1`.ice.ini.nc

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
echo =============
echo 'prepare fesom'
echo =============
setenv yeartoicemod `cat $timmsoutdir/$runid.yeartoicemod.dat`
setenv yearfromicemod `cat $timmsoutdir/$runid.yearfromicemod.dat`
setenv CLOCKFILE $runid.clock
#echo 'first check correctness of' $fesomdatadir/$CLOCKFILE
#setenv VAR11 `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $1}'`
#setenv VAR12 `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
setenv OLDYEAR `awk -F";" 'NR == 1 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
#setenv NEWTIME `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $1}'`
#setenv NEWDAY `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $2}'`
#setenv NEWYEAR `awk -F";" 'NR == 2 { print $1}' $fesomdatadir/$CLOCKFILE | awk '{print $3}'`
#echo "$VAR11"
#echo "$VAR12"
echo "$OLDYEAR"
#echo "$NEWTIME"
#echo "$NEWDAY"
#echo "$NEWYEAR"
#if ($NEWDAY != 1 || $OLDYEAR == $NEWYEAR || $OLDYEAR != $yeartoicemod || $NEWYEAR != $yearfromicemod) then
# echo "--------------------------------------------"
# echo "ERROR: invalid clock file - Mr. Timms terminates"
# echo "--------------------------------------------"
# exit
#else 
# echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# echo "$CLOCKFILE is ok: We are GO for FESOM launch sequence."
# echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#endif
echo "archive fesom output using year.yearfrac format ... "
set year = `echo $yeartoicemod | cut -d "." -f1`
set yearfrac = `echo $yeartoicemod | cut -d "." -f2`
set names = ('forcing.diag' 'ice.diag' 'ice.mean' 'mesh.diag' 'oce.diag' 'oce.mean' 'oce' 'ice')
foreach name ($names)
   mv $fesomdatadir/$runid.$year.$name.nc $fesomdatadir/$runid.$year.$yearfrac.$name.nc
   echo "$name is done"
end
set names = ('ice.ini' 'oce.ini')
foreach name ($names)
   cp $fesomdatadir/$runid.$year.$name.nc $fesomdatadir/$runid.$year.$yearfrac.$name.nc
   echo "$name is done"
end
sleep 5
setenv newmeshdir $fesommeshdir/$yearfromicemod
setenv configfile $fesomsourcedir/namelist.config
echo 'now modify namelist.config using' $OLDYEAR'(I know, its strange) and' $newmeshdir
cp $configfile $fesomrundir/namelist.config.$yeartoicemod
sed -i "s~yearnew=.*~yearnew=$OLDYEAR~" $configfile
sed -i "s~MeshPath=.*~MeshPath='$newmeshdir/'~" $configfile
sed -i "s~runid=.*~runid='$runid'~" $configfile

setenv oldmeshdir $fesommeshdir/$yeartoicemod
echo 'remove mesh/dist folder from last year'
rm -r $oldmeshdir/dist

if ($2 == 'only') then
 echo "==========================================="
 echo "exit after prepare_fesom has been completed"
 echo "==========================================="
 exit
endif

launchfesom:
echo ==========================================
echo 'now launch fesom'
echo ==========================================
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
echo ==========================================
echo 'now launch sleeping lookup job'
echo ==========================================
setenv yearfromicemod `cat $timmsoutdir/$runid.yearfromicemod.dat`
cd $masterdir
./check4fesomdata.csh $fesomdatadir $runid $yearfromicemod yes $timmsoutdir > $timmsoutdir/check4fesomdata.log &
sleep 2


echo =================================================
echo 'Mr. Timms script completed; FESOM now runs for year'
cat $timmsoutdir/$runid.yearfromicemod.dat
echo =================================================

exit

