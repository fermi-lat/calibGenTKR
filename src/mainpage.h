/** @mainpage package calibGenTKR
  @author Leon Rochester

  @section intro Introduction
  This package is set up to run the root macro doBadStripsCalib as a compiled
  program

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
  1) path to input digi.root file
            The current test file is found at: 
            glast_03/EM2003/rootFiles/em_v1r030302p5/digi
  2) name of digi.root file
            The current test file is named 
            ebf031006235353_digi.root
  3) prefix for output files
  4) Maximum number of events read

  The default input file is named options.txt and is expected 
            in the /src/test directory of the package
  An alternate filename may be passed as an argument
  @endverbatim

*/

