#include "BadStripsCalib.h"
#include "TSystem.h"
#include <iostream>
#include <string>

/* to do:
  figure out better local average, for now we use simple average
  pass string for output files
  clean out debris
  input params from text file?
  put ToT stuff on a switch; most likely needs to be moved to a different
     calib process.
  */

/** @file generateBadStrips
@brief Driver for bad strips calibration
*/
int main(int argn, char** argc) {
    
#ifdef WIN32
    gSystem->Load("libTree.dll");   
#endif

    unsigned int numEvents = 5000000;
	enum Sample {EM, TEST};
	Sample current = EM;
    bool attended = false;

    // To do: these strings will need to come from an input parameter file
    // current example is from glast_03/EM2003/rootFiles/em_v1r030302p5/digi 
    std::string sourceDirectory("c:/Glast/files/em/");
    std::string sourceFile     ("ebf031006235353_digi.root");
    std::string outputString   ("ebf031006235353");
    outputString = "example";
    std::string path;
    if( current==TEST)    { path = ::getenv("BADSTRIPSCALIBROOT");}
    else if (current==EM) { path = sourceDirectory; }

    std::string digiFileName(path);

    if (current==TEST) {
        digiFileName += "/src/test/digi.root";
    }	
    else if (current==EM) {
		digiFileName += sourceFile; 
	}
    std::cout << "digi File: " << digiFileName << std::endl;

    // not needed for this calibration
    std::string mcFileName = "";
    std::string reconFileName = "";

    if (argn > 1) mcFileName = argc[1];
    if (argn > 2) digiFileName = argc[2];
    if (argn > 3) reconFileName = argc[3];

    std::string outputPrefix(::getenv("CALIBGENTKRROOT"));
    outputPrefix += "/output/"+outputString;
    char* c_prefix = const_cast<char*>(outputPrefix.c_str());

    BadStripsCalib r(digiFileName.c_str(), reconFileName.c_str(), mcFileName.c_str(), c_prefix);
    r.Go(numEvents);
	r.Finish();
    r.WriteHist();
    
    if (attended) {
        char istop;
        std::cout << "Hit Enter key to exit: ";
        std::cin.get(istop);
    }
    return 0;
}







