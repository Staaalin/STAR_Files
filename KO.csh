#!/bin/csh

cd ~/STAR_Files
csh Delate.csh

condor_rm svianping

cd ~/STAR_Files/KFParticle4Lambda
./sub.csh


