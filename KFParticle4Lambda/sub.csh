#!/bin/csh

rm -rf ~/STAR_Files/KFParticle4Lambda/xml
mkdir ~/STAR_Files/KFParticle4Lambda/xml

set MainDir=`pwd`

set XmlDir=./xml
if(! -e $XmlDir) exit
cd $XmlDir
set SubXml=sub.xml
if(-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=2 #40
set nFileTotal=all

# print xml file
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
echo \<command\>$MainDir/run\.csh\</command\> >> $SubXml
# echo \<stdout URL=\"file:$MainDir/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
echo \<stdout URL=\"file:/star/data01/pwg/svianping/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P16id,filetype=daq_reco_PicoDst,trgsetupname~production_pAu200_2015,runnumber\[\]16124017-16159024,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P17id,filetype=daq_reco_PicoDst,trgsetupname~dAu200_production_2016,runnumber\[\]17132063-17141003,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P17id,filetype=daq_reco_PicoDst,trgsetupname~dAu62_production_2016,runnumber\[\]17141041-17148003,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
echo \<input URL=\"catalog:star\.bnl\.gov\?production=P17id,filetype=daq_reco_PicoDst,trgsetupname~dAu39_production_2016,runnumber\[\]17160011-17169018,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P17id,filetype=daq_reco_PicoDst,trgsetupname~dAu20_production_2016,runnumber\[\]17149053-17160009,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P23id,filetype=daq_reco_PicoDst,trgsetupname~production_dAu200_2021,runnumber\[\]22180043-22188007,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P16id,filetype=daq_reco_PicoDst,trgsetupname~production_pp200trans_2015,runnumber\[\]16064034-16093018,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
echo \<output fromScratch=\"root_\$JOBINDEX\.log\" toURL=\"file:/star/data01/pwg/svianping/log/\" /\> >> $SubXml
echo \<output fromScratch=\"output_\$JOBINDEX\.root\" toURL=\"file:/star/data01/pwg/svianping/output/\" /\> >> $SubXml
echo \<output fromScratch=\"KFParticleQA_\$JOBINDEX\.root\" toURL=\"file:/star/data01/pwg/svianping/output/\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml
