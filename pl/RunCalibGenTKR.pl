#!/usr/local/bin/perl

use strict;

my $jobFile = "job.dat";

my $cmtPath = '/nfs/farm/g/glast/u06/chen/calib:/nfs/farm/g/glast/u06/chen/svac/:/nfs/farm/g/glast/u10/builds/EngineeringModel/EngineeringModel-v3r0402p7/';
my $cmtDir = '/nfs/farm/g/glast/u06/chen/calib/calibGenTKR/v0r1p2/cmt';
my $exeDir = '/nfs/farm/g/glast/u06/chen/calib/calibGenTKR/v0r1p2/rh9_gcc32';
my $digiRootDir = '/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/digi/';
my $outputDir = '/nfs/farm/g/glast/u06/chen/temp/';

open(JOBFILE, "<$jobFile") || die "Can't open $jobFile for input, abortted!";

while(<JOBFILE>) {

    my $outputPrefix;
    ($outputPrefix) = split;
    if( $outputPrefix eq "") { last; }

    my $digiRootFile = $outputPrefix.'_digi.root';

    my $dumpFile = $outputDir.'calibGenTKR_'.$outputPrefix.'.dump';

#create option xml file for doBadStripsCalib.exe

    my $optionFile = $outputDir.'calibGenTKR_'.$outputPrefix."_opt.dat";
    open(OPTIONFILE, ">$optionFile") || 
	die "Can't open $optionFile for input, abortted!";

    print OPTIONFILE qq{<?xml version="1.0" ?> \n};
    print OPTIONFILE qq{<!-- unix test version --> \n};

    print OPTIONFILE qq{<!DOCTYPE ifile SYSTEM "\$(CALIBGENTKRROOT)/xml/ifile.dtd" > \n};

    print OPTIONFILE qq{<ifile  cvs_Header="\$Header: /nfs/slac/g/glast/ground/cvs/calibGenTKR/pl/RunCalibGenTKR.pl,v 1.1 2004/06/18 22:32:14 xchen Exp $" cvs_Revision="\$Revision: 1.1 $"> \n};
    print OPTIONFILE qq{  <section name="parameters"> input parameters for TKR bad strips calibration \n};
    
    print OPTIONFILE qq{    <item name="detectorType" value="EM1"> valid types are EM1, EM2, LAT_2Towers and LAT_Full</item> \n};

    print OPTIONFILE qq{    <item name="maxOccupancy" value = "0.01"> occupancy above ambient to qualify as hot</item> \n};

    print OPTIONFILE qq{    <item name="sourceFilePath" value= "$digiRootDir" >  common part of the path to the digi.root files (may be blank) </item> \n};
    print OPTIONFILE qq{    <item name="sourceFileList" value= "$digiRootFile" > blank-delimited list of input files, including any part of the path not specified above </item> \n}; 

    print OPTIONFILE qq{    <item name="xmlPath" value= "$outputDir">  path to location of output xml files </item> \n};
    print OPTIONFILE qq{    <item name="histPath" value= "$outputDir">  path to location of output hist files </item> \n};
    print OPTIONFILE qq{    <item name="outputPrefix"    value= "calibGenTKR_$outputPrefix" >  prefix for output files </item> \n};
    print OPTIONFILE qq{    <item name="numEvents"       value= "500000" > max number of events to read </item> \n};
    print OPTIONFILE qq{  </section> \n};
    print OPTIONFILE qq{</ifile> \n};

    close(OPTIONFILE);

#create shell file to execute doBadStripsCalib.exe
    my $shellFile = $outputDir.'calibGenTKR_'.$outputPrefix.".scr";
    open(SHELLFILE, ">$shellFile") || 
	die "Can't open $shellFile for input, abortted!";
    print SHELLFILE qq{#!/bin/csh \n \n};
    print SHELLFILE qq{unsetenv LD_LIBRARY_PATH \n};
    print SHELLFILE qq{setenv CMTPATH $cmtPath \n};
    print SHELLFILE qq{source $cmtDir/setup.csh \n};
    print SHELLFILE qq{$exeDir/doBadStripsCalib.exe $optionFile \n};
    close(SHELLFILE);
    system("chmod +rwx $shellFile");
#    system("$shellFile");
    system("bsub -q short -o $dumpFile $shellFile");
}

close(JOBFILE);
