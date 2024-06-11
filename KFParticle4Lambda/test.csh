#!/bin/csh
#starver SL19b
# starver DEV
# starver SL20a

## AuAu@27GeV
# # inputs
# set iJob=$1
# # settings
# set nRun=18
# set mRun=Run${nRun}
# set mEnergy=27.0
# set ListDir=./datalist/
# set MainDir=`pwd`

## pAu@200GeV or dAu@200GeV
# inputs
set iJob=$1
# settings
set nRun=16
set mRun=Run${nRun}
set mEnergy=62.0
set ListDir=./datalist/
set MainDir=`pwd`

## AuAu@200GeV
# # inputs
# set iJob=$1
# # settings
# set nRun=18
# set mRun=Run${nRun}
# set mEnergy=200.0
# set ListDir=./datalist/
# set MainDir=`pwd`

echo $mEnergy $nRun

set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/$mEnergy.list.`printf "%.6d" ${iJob}`
#set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/test.list

root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$iJob,$nRun,$mEnergy,\"$ListDir\"\)
# /afs/rhic.bnl.gov/star/packages/SL23c/.sl73_gcc485/bin/root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$iJob,$nRun,$mEnergy,\"$ListDir\"\)
