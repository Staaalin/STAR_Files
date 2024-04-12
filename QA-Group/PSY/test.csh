#!/bin/csh
#starver SL19b
# starver DEV
starver SL20d

source /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh
cons

## pAu@200GeV or dAu@200GeV
# inputs
set iJob=$1
# settings
set nRun=18
set mRun=Run${nRun}
set mEnergy=26.5
set ListDir=datalist/
set MainDir=`pwd`

set cen = 1
set opt_weight = 1

echo $mEnergy $nRun

set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/$mEnergy.list.`printf "%.6d" ${iJob}`

cd /star/u/svianping/STAR_Files/QA-Group/StRoot

root4star -b -q /star/u/svianping/STAR_Files/QA-Group/PSY/Gamma_QA.C\($cen,$opt_weight,\"/star/u/svianping/STAR_Files/QA-Group/PSY/$FILELIST\"\)
