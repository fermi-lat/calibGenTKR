/** @mainpage package calibGenTKR
  @author Leon Rochester

  @section intro Introduction
  This package is set up to run the root macro doBadStripsCalib as a compiled
  program.

  There are several running modes. The current setting is simple cuts with no
  count check. The other modes need a bit of work.

  ToDo:
  <br>Use TChain for input, even if only one file

  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

  @section input Structure of input file
  @verbatim
  File consists of four lines:
  1) Path to input digi.root file
            The current test file is found at: 
            /nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/digi/
  2) Name of digi.root file
            The current test file is named: 
            ebf031006235353_digi.root
  3) Prefix for output files. Appended to this is:
            _deadStrips.xml, _hotStrips.xml or _hist.root
  4) Maximum number of events to be processed
            Set to a large number to get all the events.
            For an EM run, about 100K events are needed.
            For the flight instrument, about 50K events/tower are needed.

  The default input file is named options.txt and is expected 
            in the /src/test directory of the package
  An alternate filename may be passed as an argument
  @endverbatim

*/

