#! /opt/star/bin/perl -w

use File::Basename;
use Getopt::Std;
use Cwd 'abs_path';     # aka realpath()
my $resubmit = "resub.sh";
my $remove = "remove.sh";
my $chopCondor = "Sample";
@jobDirs = (
#"/star/data01/pwg/zhiwanxu/ESE-11GeV/P23ia_EPD1_run2/cen0"
#"/star/data01/pwg/zhiwanxu/ESE-11GeV/P23ia_EPD1/cen0"
#"/star/data01/pwg/zhiwanxu/QA-group/11GeV/P23ia_QA/cen0"
"/star/data01/pwg/zhiwanxu/QA-group/11GeV/P23ia_cen/cen0"
#"/star/data01/pwg/zhiwanxu/ESE-11GeV/P23ia_EPD_v2/cen0"
);

@flowCents = (
"cen0",
#"cen1","cen2","cen3","cen4","cen5","cen6","cen7","cen8","cen9"
);

foreach $jobDir (@jobDirs) {print `rm -rf $jobDir/*remove* \n`;}

foreach $eachCent (@flowCents) {

$totalJobs=0.;
$finishedJobs=0.;
    foreach $jobDir (@jobDirs) {
       if (-e "$jobDir/$eachCent.$remove") { print `rm $jobDir/$eachCent.$remove \n`;}
       if (-e "$jobDir/$eachCent.$resubmit") { print `rm $jobDir/$eachCent.$resubmit \n`;}
       if (-e "$jobDir/$eachCent.$chopCondor") { print `rm -rf $jobDir/$eachCent.$chopCondor \n`;}
       print `mkdir $jobDir/$eachCent.$chopCondor \n`;

   foreach (glob("$jobDir/sched*_*_*.condor")){
	my $eachScript = $_;
        open (prototypeMacro, "$eachScript");
       chomp $eachScript;
       @fields = split(/\//,$eachScript) ;
       $eachScriptNoPath= $fields[$#fields];
       ($eachJOBID = $eachScriptNoPath) =~ s/\.condor// ;
       chomp $eachJOBID;
       @names = split(/_/,$eachJOBID) ;
       $eachNumber= $names[$#names];
#         print "$eachNumber\n";
         $eachName = $names[0];
#         print "$eachName\n";
         
    $count =0.;
    $number = 13;
    while ($eachLine = <prototypeMacro>) {
      
          if(($count % $number)==0) {open (macroFile,">$jobDir/$eachCent.$chopCondor/$eachName\_$eachNumber.condor");}
         
         print macroFile $eachLine;
         $count = $count+1;
         if(($count % $number)==0) {close macroFile; $eachNumber = $eachNumber -1;}
          
          
     }
     
     close prototypeMacro;

}

print  `ls $jobDir/$eachCent.$chopCondor/ |wc -l`;
print "condor files have been chopped. Remeber to examine if the lines match!\n";

foreach (glob ("$jobDir/sched*.csh") ){
      my $eachScript = $_;
      chomp $eachScript;
      @fields = split(/\//,$eachScript) ;
      $eachScriptNoPath= $fields[$#fields];
      ($eachJOBID = $eachScriptNoPath) =~ s/\.csh// ;
      chomp $eachJOBID;


#if ( (!(-s "$jobDir/$eachJOBID.$eachCent.gamma112_nop_EP11_Boost0.root")) ) {
# if ( (!(-s "$jobDir/$eachJOBID.$eachCent.gamma112_nop_EP1_Boost0.root")) ) {
#if ( (!(-s "$jobDir/$eachJOBID.$eachCent.gamma112_pipi_EP0_Boost0.root")) ) {
if ( (!(-s "$jobDir/$eachJOBID.$eachCent.centrality.def.root")) ) {
#if ( (!(-s "$jobDir/$eachJOBID.$eachCent.v2.root")) ) {
##if ( (!(-s "$jobDir/$eachJOBID.$eachCent.gamma112_EPD2_EP11_Boost0.root")) ) {
#	print "$eachJOBID\n";
#    	if (-e "$jobDir/$eachCent/$eachJOBID.run.log") { print `rm $jobDir/$eachCent/$eachJOBID.run.log \n`; }
	$removecomd = "rm $eachJOBID.*.root";
	$submitCommand = "cd $jobDir;cp $jobDir\/$eachCent.$chopCondor\/$eachJOBID.condor .; condor_submit $jobDir\/$eachJOBID.condor";
#     	$submitCommand = `grep "condor_submit" $eachScript`;
#     	chomp $submitCommand;
#     	$submitCommand =~ s/\#//;

#	if ($submitCommand=~/direct/) { #if at RCF, /direct/star+dataxx/ does not work ! ANNOYING
#		$submitCommand=~ s/\/direct//;
#		$submitCommand=~ s/\+/\//;   # change /direct/star+dataxx to /star/dataxx
#      	}

#	$submitCommand=~ s/\]\" /\]\" < /;   # change /direct/star+dataxx to /star/dataxx

	open(macroFile,">>$jobDir/$eachCent.$resubmit");
	print macroFile "$submitCommand \n";

        open(macroFile2,">>$jobDir/$eachCent.$remove");
        print macroFile2 "$removecomd \n";

  }

}


   my @tempScripts = glob( "$jobDir/sched*.csh" );
 my @tempFinished = glob( "$jobDir/sched*$eachCent.centrality.def.root" ) ;
#my @tempFinished = glob( "$jobDir/sched*$eachCent.v2.root" ) ;

              $totalJobs +=scalar(@tempScripts);
	   $finishedJobs +=scalar(@tempFinished);

    }

my $perCentFinished = 0.;
if ($totalJobs != 0) {$perCentFinished = ($finishedJobs/$totalJobs)*100.;}
 print "for $eachCent, finished $finishedJobs jobs out of $totalJobs, ".$perCentFinished."% completed \n";
}
exit;

