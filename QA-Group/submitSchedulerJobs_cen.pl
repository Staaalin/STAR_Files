#! /opt/star/bin/perl -w

use File::Basename;
use Getopt::Std;
use Cwd 'abs_path';     # aka realpath()

my %opt; 
getopts('htoz:b:m:adn:q:x:y:s',\%opt);
my $usage = "$0  configfile combinedDenominatorDir\n";
$usage .= "\t if no combinedDenominatorDir, assuming it is the first time submition \n";
#----

my $configFile = shift or die $usage;
my $combinedDenomDir = shift;

$configFile =~ /.config/ or die "are you sure that $configFile is a config file ? \n";

my $pwd = `pwd`;
if ($pwd=~/direct/) { #if at RCF, /direct/star+dataxx/ does not work ! ANNOYING
	$pwd=~ s/\/direct//;
	$pwd=~ s/\+/\//;   # change /direct/star+dataxx to /star/dataxx
}
chomp($pwd);

# read in config file
  open (CONFIG, "< $configFile") or die "can't open $configFile $!";
  while (<CONFIG>){
       chomp;                          # no newline
       s/#.*//;                        # no comments
       s/^\s+//;                       # no leading white
       s/\s+$//;                       # no trailing white
       next unless length;             # anything left?
       my ($var, $value) = split (/\s*:\s*/, $_,2);
       $ConfigMap{$var} = $value;
  }
  close(CONFIG);

@flowCents = (
 "cen0"
);
@RunCents =( 
"cen0"
#"cen1", "cen2", "cen3", "cen4", "cen5", "cen6", "cen7", "cen8", "cen9"
);

my $macroName        = $ConfigMap{"macroName"};
my $codeName	     = $ConfigMap{"codeName"};
my $trimPerl	     = $ConfigMap{"trimPerl"};
my $libVersion       = $ConfigMap{"libVersion"};
my $runJob           = $ConfigMap{"scriptName"};
my $jobDir           = $ConfigMap{"jobDir"};
my $tmpDir	     = $ConfigMap{"tmpDir"};
my $jobName          = $ConfigMap{"jobName"};
my $logFile          = $ConfigMap{"logFileName"};
my $cumulFile        = $ConfigMap{"histoFile"};
my $schedName        = $ConfigMap{"schedName"};
my $maxFilesPerJob   = $ConfigMap{"maxFiles"};

# move the previouse run as old
if ($combinedDenomDir) {
    if (-e "$jobDir/$jobName"){
 print `mv $jobDir/$jobName   $jobDir/$jobName."old"`;
    }
}

# make job directories
 print `mkdir $jobDir/$jobName`;

# copying source code to $jobName
print `cp -Lr $jobDir/StRoot $jobDir/$jobName`;

print `cp $configFile $jobDir/$jobName/$configFile`;

# cons
 my $consFileName="consDir";
 open ( consScript, ">".$consFileName.".csh");
 print consScript  "#!/bin/csh \n";
 print consScript  "cd  $jobDir/$jobName \n";
 print consScript  "source \$GROUP_DIR/.starver $libVersion \n";
# print consScript  "stardev \n";
# print consScript  "cons \n";
 close consScript;

 print "$consFileName.csh has been generated \n";
 print "now run $consFileName.csh  : \n";
 print `csh -x ./$consFileName.csh`;
 print `rm ./$consFileName.csh`;


# making centrality directories, and lib links
 foreach $cent(@flowCents) {
 if(-e "$jobDir/$jobName/$cent") { print `rm -r $jobDir/$jobName/$cent \n`;}
 print `mkdir $jobDir/$jobName/$cent \n`;
 if (-e "$jobDir/$codeName"){
 print `cp $jobDir/$codeName  $jobDir/$jobName/$cent \n`;
 }
 if (-e "$jobDir/$macroName"){
 print `cp $jobDir/$macroName $jobDir/$jobName/$cent \n`;
 }
 if (-e "$jobDir/$trimPerl"){
 print `cp $jobDir/$trimPerl  $jobDir/$jobName/$cent \n`;
 }
 print `ln -s $jobDir/.sl73_gcc485 $jobDir/$jobName/$cent/.sl73_gcc485 \n`;

# producing shell script

 open (shellScript,">$jobDir/$jobName/$cent/$runJob");
  print shellScript "#! /bin/csh \n";
  print shellScript "set logfile=$logFile \n";
  print shellScript "set prototypeWDIR=$jobDir/$jobName/$cent \n";
  print shellScript "\n";
  print shellScript "set WDIR=$tmpDir/\$JOBID \n"; 
  print shellScript " \n";
  print shellScript "if ( -e \"$WDIR\" ) then \n";
  print shellScript "rm -rf $WDIR \n";
  print shellScript "endif \n";
  print shellScript " \n";

  print shellScript "mkdir -p \$WDIR || exit 1 \n";
  print shellScript "cd \$WDIR \n";

  print shellScript "if ( -e \"\$prototypeWDIR/$codeName\" ) then \n";
  print shellScript "cp -Lr \$prototypeWDIR/$codeName . \n";
  print shellScript "endif \n";

  print shellScript "if ( -e \"\$prototypeWDIR/$macroName\" ) then \n";
  print shellScript "cp -Lr \$prototypeWDIR/$macroName . \n";
  print shellScript "endif \n";

  print shellScript "if ( -e \"\$prototypeWDIR/$trimPerl\" ) then \n";
  print shellScript "cp -Lr \$prototypeWDIR/$trimPerl . \n";
  print shellScript "endif \n";

  print shellScript "cp -Lr \$prototypeWDIR/.sl73_gcc485 . \n";

  print shellScript "cp \$FILELIST tes.list \n";
  print shellScript "perl $trimPerl \n";
  print shellScript "mv test1.list test.list \n";

  print shellScript  "setenv NODEBUG yes \n";
  print shellScript  "source \$GROUP_DIR/.starver $libVersion \n";

  foreach $RunCent(@RunCents) {
        @fields = split(/cen/,$RunCent) ;
        $centrality= $fields[$#fields];
        print `echo $centrality \n`;

	print shellScript "root -b -q '$macroName++' >& \${logfile} \n";
	print shellScript "if ( -e \"cen$centrality.$cumulFile\") then \n";
  	print shellScript "mv -f cen$centrality.$cumulFile  \$prototypeWDIR/sched\$JOBID.cen$centrality.$cumulFile \n";
  	print shellScript "gzip $logFile \n";
  	print shellScript "mv -f $logFile.gz    \$prototypeWDIR/sched\$JOBID.$logFile.gz \n";
  	print shellScript "else \n";
  	print shellScript "mv -f $logFile    \$prototypeWDIR/sched\$JOBID.cen$centrality.$logFile \n";
  	print shellScript "echo \"Action did not produce the expected $cumulFile file\" \n";
  	print shellScript "endif \n";
}
  print shellScript "rm -fr \$WDIR \n";
  print shellScript "exit 0 \n";
 close shellScript;

 }

my $xmlFile;

  foreach $cent(@flowCents) { 
     open (macroFile,">$jobDir/$jobName/$cent/$schedName");
     print macroFile "<?xml version=\"1.0\" encoding=\"utf-8\" ?> \n";
     print macroFile "<job simulateSubmission =\"false\" maxFilesPerProcess =\"$maxFilesPerJob\" fileListSyntax=\"xrootd\"> \n";
     print macroFile "<command>$jobDir/$jobName/$cent/$runJob</command> \n";
     print macroFile "<stdout URL=\"file:$jobDir/$jobName/$cent/\$JOBID.out\" /> \n";
#	print macroFile "<input URL=\"catalog:star.bnl.gov?trgsetupname=production_11p5GeV_2020,production=P23ia,filetype=daq_reco_picoDst,filename~st_physics,runnumber[]21020007-21041000,storage=local\" nFiles=\"all\" /> \n";
	print macroFile "<input URL=\"catalog:star.bnl.gov?trgsetupname=production_11p5GeV_2020,production=P23ia,filetype=daq_reco_picoDst,filename~st_physics,storage=local\" nFiles=\"all\" /> \n";
	#print macroFile "<input URL=\"catalog:star.bnl.gov?trgsetupname=production_11p5GeV_2020,production=P23ia,filetype=daq_reco_picoDst,filename~st_physics,storage=nfs\" nFiles=\"all\" /> \n";
     print macroFile "</job> \n"; 

   close macroFile;
 }

 my $SubmitJobName="submitSched.csh";

 foreach $cent(@flowCents) {

    open ( submitScript, ">$jobDir/$jobName/$cent/$SubmitJobName");
     print submitScript "#! /bin/csh \n";
     print submitScript "cd $jobDir/$jobName/$cent \n";
     print submitScript "star-submit $schedName \n";
     print submitScript "exit 0 \n";
    close submitScript;
 }

 foreach $cent(@flowCents) { 
 	print `chmod +x $jobDir/$jobName/$cent/$schedName \n`;
 	print `chmod +x $jobDir/$jobName/$cent/$runJob \n`;
 	print `chmod +x $jobDir/$jobName/$cent/$SubmitJobName \n`;
 	print `csh -x $jobDir/$jobName/$cent/$SubmitJobName \n`;
 }

exit;




#perl parallFlowJobSubmit.pl --d /auto/pdsfdv36/starebye/aihong/MinbiasP01hi/copiedFromPdsfdv08 --i 5 --o test


####$pwd = abs_path ( $ENV { 'PWD' } );
