/** @mainpage package calibGenTKR
  @author Leon Rochester

  @section intro Introduction
  This package is set up to run the root macro doBadStripsCalib as a compiled
  program.

  @section operation Operation

  A large number of events are read in from one or more digi.root files. The events should
  be reasonably populated with real tracks. Otherwise, the process becomes sensitive to 
  the presence of noise hits, and the results will start to depend on the details of the 
  threshold settings. (If this becomes a problem in the full LAT, we can do something fancier
  to handle it. Noise hits are useful to detect hot strips, but complicate the analysis for 
  dead strips.) 
  
  In principle, each strip should be struck by some minimum number of tracks,
  perhaps 40. For a single tower, with cosmic trigger, this implies about 50K events;
  for the full LAT, 800K.
  The hit strips are accumulated in histograms, one for each plane.

  In addition, any hits in illegal towers or planes are accumulated into two histograms whose
  contents should be zero.

  The histograms are analyzed in one of several ways, only the first of which 
  has been thoroughly tested:

  <ul> 
  <li>Simple: the strips with zero hits are flagged as dead; those with occupancy greater than 
  <em>maxOccupancy</em> above the ambient (default is 1%)
  are flagged as hot. For events in a full-sized tower, this scheme should be fine.</li>

  <li>Standard: (This option needs some work.) The number of hits for each 
  strip is tested against the expected value,
  based on the average number of hits/strip for that plane. Strips with counts more than 
  <em>nSigma</em> (default
  is 6) sigma
  below expectation, or zero counts, are flagged as dead; those with occupancy greater than
  <em>hotFactor</em> (default is 10) times the ambient are flagged as hot.</li>

  <li>Fit to local average: (This one needs even more work!) Like the previous, 
  but the average used is the value of a polynomial fit
  to the data in the histogram. This is supposed to compensate for effects in the EM, where
  there are many fewer counts at the edges of the distributions.</li>
  </ul>

  A further test can be made on the dead strips, as to whether the expected value, based on the 
  appropriate average, would have been high enough to preclude accidental zeroes. 
  (This test is turned off by default.)

  The job produces two xml files, one containing the dead strips, 
  and another containing the hot strips.
  In addition, the histograms are output in a root file.

  @section input Input parameters

  The input parameters are supplied in an xml file.
  The default example file is options.xml in the /src/test directory, and accesses
  a unix input file of EM1 data. An alternate file may be passed as an argument to 
  the program.
  
  @param detectorType
  Sets number of layers, towers, tower numbers. Allowed values are:
  EM1, EM2, LAT_2Towers, LAT_Full.
  @param maxOccupancy
  Sets threshold for calling a strip hot. Default is 0.01 above ambient
  @param sourceFilePath
  The common part of the paths to the input digi.root files. Can be set to "",
  or omitted, if each input file is in a different place.
  @param sourceFileList
  Blank-delimited list of input files, including that part of the path not included
  in sourceFilePath. The example file is found at
  /nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/digi/.
  @param xmlPath
  the location of the output xml files. Default is the output directory of calibGenTKR
  @param histPath
  the location of the output histogram file. Default is the output directory of calibGenTKR
  @param outputPrefix
  Prefix for the output files. This string is prepended to the names:
  _deadStrips.xml, _hotStrips.xml, and _hist.root.
  @param numEvents
  Maximum number of events to be processed. The actual number is the 
  minimum of this number and the number of events in the TChain.

  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements



*/
