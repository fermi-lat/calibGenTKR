/** @mainpage package calibGenTKR
  @author Leon Rochester

  @section intro Introduction
  This package is set up to run the root macro doBadStripsCalib as a compiled
  program.

  There are several running modes. The current setting is simple cuts with no
  count check. The other modes need a bit of work.

  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

  @section input Input parameters

  The input parameters are supplied in an xml file.
  The default file is options.xml in the /src/test directory, and accesses
  a unix input file. An alternate file may be passes as an argument to 
  the test program.
  
  @param sourceFilePath
  Path to the input digi.root files
  @param sourceFileList
  Blank-delimited list of input files. The test file is found at
  /nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/digi/.
  @param outputPrefix
  Prefix for the output file. This string is prepended to the names:
  _deadStrips.xml, _hotStrips.xml, and _hist.root.
  @param numEvents
  Maximum number of events to be processed. The actual number is the 
  minimum of this number and the number of events in the TChain.
*/
