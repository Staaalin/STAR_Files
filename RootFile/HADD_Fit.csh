#!/bin/csh

set InputName = "~/Result/HADDrR_0.root"
set OutputName = "~/Result/cHADDrR_0.root"

root4star -q -b HADDr_Fit.C\($InputName,$OutputName\)