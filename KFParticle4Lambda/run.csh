#!/bin/csh

# AuAu RUN=18 Energy=27.0 

# settings
set nRun=16
set mRun=Run${nRun}
set mEnergy=200
set ListDir=/star/u/svianping/STAR_Files/KFParticle4Lambda/datalist #TODO
set MainDir=/star/u/svianping/STAR_Files/KFParticle4Lambda #TODO
set TempDir=/home/tmp/svianping #TODO
# inputs
# set JOBINDEX=$1
# set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/$mEnergy.list.`printf "%.6d" ${JOBINDEX}`
#set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/test.list
#set JOBID=ScriptTestSandbox

echo $FILELIST

set WorkDir=${TempDir}/$JOBID
mkdir -p $WorkDir
cd $WorkDir
cp -Lr $MainDir/setDEV2.csh .
cp -Lr $MainDir/lMuDst.C .
source setDEV2.csh
cp -Lr $MainDir/readPicoDst.C .
cp -Lr $MainDir/.sl73_x8664_gcc485 .

set RootLog=root_${JOBINDEX}.log
if(-e $RootLog) rm $RootLog

root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$JOBINDEX,$nRun,$mEnergy,\"$ListDir\"\) >& $RootLog

set Iter=0
while( `grep -sc '(ret%10)<=kStFatal' $RootLog` )
	@ Iter++
	if( $Iter > 5 ) break
	rm $RootLog
	rm *.root
	root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$JOBINDEX,$nRun,$mEnergy,\"$ListDir\"\) >& $RootLog
end

# mv *.log  $MainDir/log/.
# mv *.log  /star/data01/pwg/svianping/log/.
# mv *.root $MainDir/output/.
# mv *.root /star/data01/pwg/svianping/output/.