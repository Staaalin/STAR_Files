#!/bin/csh
#starver SL19b
# starver DEV
rm -r .sl73_gcc485
setenv NODEBUG yes
starver SL20d

echo "If Source? 0:[NO] 1:[YES]"
set IfSource = "$<"
if ($IfSource == 1) then
    cp /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh /star/u/svianping/STAR_Files/QA-Group/PSY/
    source setDEV2.csh
endif
echo "If Cons? 0:[NO] 1:[YES]"
set IfCons = "$<"
if ($IfCons == 1) then
    cons
endif


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

cd /star/u/svianping/STAR_Files/QA-Group/PSY

# root4star -b -q /star/u/svianping/STAR_Files/QA-Group/PSY/Gamma_QA.C\($cen,$opt_weight,\"/star/u/svianping/STAR_Files/QA-Group/PSY/$FILELIST\"\)
root4star -b -q /star/u/svianping/STAR_Files/QA-Group/PSY/RunAnalyzer_QA.C\($cen,$opt_weight,\"/star/u/svianping/STAR_Files/QA-Group/PSY/$FILELIST\"\)
