#!/bin/bash

cd output
mkdir -p temp
cp ../hadd2.pl ./temp/

cp KFParticleQA_* ./temp/
cd temp
perl hadd2.pl
rm KFParticleQA_*
perl hadd2.pl
