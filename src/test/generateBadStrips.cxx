#include "badStripsCalib.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>

/* to do:
  figure out better local average, for now we use simple average
  pass string for output files
  input params from text file?
  */

/** @file generateBadStrips
@brief Driver for bad strips calibration
*/

/* attempt to do a bit of processing on the options file
// doesn't work, don't know why
std::string getString(std::string& inString) {
    std::string temp = inString;
    int npos;
    npos = temp.find("X");
    if(npos>-1) temp = temp.substr(0, npos-1);
    npos = temp.find_first_not_of(" ");
    if(npos>-1) {
        temp = temp.substr(npos);
    } else {
        temp = "";
        return temp;
    }

    npos = temp.find_last_not_of(" ");
    if(npos>-1) temp = temp.substr(0,npos+1);
    return temp;
}

std::string nextString(std::ifstream& file) {
    std::string temp("");

    while (temp=="") {
        std::getline(file, temp);
        std::cout << "*" << temp << "*" << std::endl;
        temp = getString(temp);
    }
    return temp;
}
*/

int main(int argn, char** argc) {
    
#ifdef WIN32
    gSystem->Load("libTree.dll");   
#endif

    bool attended = false;
    // current example is from glast_03/EM2003/rootFiles/em_v1r030302p5/digi/ 
    std::string sourceDirectory("c:/Glast/files/em/");
    std::string sourceFile     ("ebf031006235353_digi.root");
    std::string outputString   ("ebf031006235353");
    std::string path = ::getenv("CALIBGENTKRROOT");
    unsigned int numEvents = 5000000;


    std::string digiFileName   (path+sourceFile);
    std::string outputPrefix(path);
    outputPrefix += "/output/"+outputString;

    std::ifstream inputFile;
    std::string tempString;
    if(argn > 1) {
        inputFile.open(argc[1], ios::in);
    }
    else {
        inputFile.open((path + "/src/test/options.txt").c_str(), ios::in);
    }

    inputFile >> sourceDirectory;
    std::cout << "Source Directory: " << sourceDirectory << std::endl;
    inputFile >> sourceFile;
    std::cout << "Source File: " << sourceFile << std::endl;

    digiFileName = sourceDirectory+sourceFile;

    inputFile >> outputString;
    std::cout << "Output string: " << outputString << std::endl;

    inputFile >> numEvents;
    std::cout << "Maximum number of events to process: " << numEvents << std::endl;

    // not needed for this calibration
    std::string mcFileName = "";
    std::string reconFileName = "";

    outputPrefix = path+"/output/"+outputString;
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







