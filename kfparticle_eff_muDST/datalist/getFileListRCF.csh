#!/bin/csh

#input begin
set energy=19.6
set nRun=19

#run18List27
set outList=run${nRun}List${energy}.list
if(-e $outList)rm -r $outList

get_file_list.pl -keys path,filename -delim '/' -cond "production=P21ic,filetype=daq_reco_picoDst,trgsetupname~production_19GeV_2019,sanity=1,tpx=1,storage!=hpss,filename~st_physics" -limit 0 > $outList
sed -i 's/^/root:\/\/xrdstar.rcf.bnl.gov:1095\//' $outList