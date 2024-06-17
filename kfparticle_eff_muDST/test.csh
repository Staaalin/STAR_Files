#!/bin/csh
#starver SL19b

# inputs
set iJob=$1
# settings
set nRun=18
set mRun=Run${nRun}
set mEnergy=27
set ListDir=./datalist/
set MainDir=`pwd`

echo $mEnergy $nRun
echo $iJob

#set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/$mEnergy.list.`printf "%.6d" ${iJob}`
set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/test.list

root4star -b -q ./analysis.C\(\"$FILELIST\",$iJob\)
