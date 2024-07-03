setenv NODEBUG yes

starver .DEV2
source $STAR/setupDEV2.csh
setenv NODEBUG yes
starver TFG21g
#setup gcc 4.4.7
setup 64b
setenv STARFPE NO

# revert to Root 6
setup 64b
setup root 6.16.00

