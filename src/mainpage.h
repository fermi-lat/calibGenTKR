/** @mainpage package calibGenTKR
  @author Leon Rochester, Xin Chen

  @section intro Introduction
  This package is set up to run two calibration applications:
  @li the root macro doBadStripsCalib as a compiled program to calibrate hot/dead strips in TKR.
  @li TOT calibration using cosmic ray muon data

  @section doBadStripCalib doBadStripCalib

  @subsection operation Operation

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

  @subsection input Input parameters

  The input parameters are supplied in an xml file.
  The default example file is options.xml in the /src/test directory, and accesses
  a unix input file of EM1 data. An alternate file may be passed as an argument to 
  the program.
  
  @param detectorType
  Sets number of layers, towers, tower numbers. Allowed values are:
  EM1, EM2, LAT_2Towers, LAT_Full. Using nPlanes and towerNumbers, below, allows more flexibility.

  @param nPlanes
  Sets the number of planes in the configuration. Default is 8. When present in the options file,
  this number overrides the detectorType.

  @param towerNumbers
  A vector of integers which specifies which towers are present in a setup. Default is "{0}",
  meaning only tower 0 is present. For the 2-tower LAT we would need "{8,9}". When present in the options file,
  this number overrides the detectorType.

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

  @section doMuonCalibTot doMuonCalibTot

  @subsection Purpose Purpose

  Calibrate potential gain differences among front end cards (or individual strips) using cosmic ray muon data.

  @subsection operation Operation

  The calibration procedure is as the following:
  @li A large number of events are read in from one or more digi.root files and corresponding recon.root files. 
  @li The application calibrates TOT using constants obtained from charge injection, the TOT value is further corrected using reconstructed muon direction, then it saves corrected TOT values in a histogram for each front end card.
  @li The histogram is fit with a Landau function whose peak is used to determine the gain value in each front end card

  @subsection jobOption Job Option file

  The application takes its first argument as the job option file. Otherwise, it uses ../src/muonCalibTot/totCalibChain_option.dat as the default.

  @li First line in the job option file is a list of digi root files separated by blank spaces.
  @li Second line in the job option file is a list of recon root files separated by blank spaces.
  @li Third line is name of output ascii file containing fitted results of the Landau distribution for each front end card in each row. First column of each row indicates whether the value is corrected using charge injection data. Second column is the plane number. Third column is the view number. Fourth column is the front end card number. Fifth column is the fitted peak value. Sixth column is the error on the peak value. Seventh column is the fitted width value. Eighth column is the error on the width value.
  @li Fourth line is name of output root file which contains one TOT histogram for every front end card.

  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements



*/
