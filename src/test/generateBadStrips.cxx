#include "badStripsCalib.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>

/** @file generateBadStrips.cxx
@brief Driver for bad strips calibration
*/

int main(int argn, char** argc) {
    
    bool attended = false;
    std::string temp;
 
    std::string sourceFilePath;
    std::string sourceFile;
    std::string path = ::getenv("CALIBGENTKRROOT");
    unsigned int numEvents = 5000000;

    std::string xmlPath(path+"/output/");
    std::string histPath(path+"/output/");
    std::string outputPrefix("test");
 
    std::string xmlFile(path+"/src/test/win_options.xml");

    if(argn > 1) {
        xmlFile = argc[1];
        std::cout << "Reading in user-specified options file: " <<  xmlFile 
            << std::endl << std::endl;
    }
    
    xml::IFile myFile(xmlFile.c_str());

    if (myFile.contains("parameters","sourceFilePath")) {
        temp = myFile.getString("parameters", "sourceFilePath");
        sourceFilePath = stripBlanks(temp);
    }

    std::cout << "Sourcefile path: " << sourceFilePath << std::endl;

    std::string sourceFileString;
    if (myFile.contains("parameters","sourceFileList")) {
        sourceFileString = myFile.getString("parameters", "sourceFileList");
    }

    std::vector <std::string> token;
    facilities::Util::stringTokenize(sourceFileString, ";, ", token);
    unsigned int nFiles = token.size();
    TChain* digiChain = new TChain("Digi");

    unsigned i;
    std::cout << "Input files:" << std::endl;
    for (i=0; i<nFiles; ++i) {
        if (token[i]=="") break;
        digiChain->Add((sourceFilePath+token[i]).c_str());
        std::cout << "   " << i+1 << ") " << token[i] << std::endl;
    }
   
    if (myFile.contains("parameters","xmlPath")) {
        temp = myFile.getString("parameters", "xmlPath");
        if (temp!="") xmlPath = temp;
    }

    if (myFile.contains("parameters","histPath")) {
        temp = myFile.getString("parameters", "histPath");
        if (temp!="") histPath = temp;
    }
    if (myFile.contains("parameters","outputPrefix")) {
        temp = myFile.getString("parameters", "outputPrefix");
        if (temp!="") outputPrefix = temp;
    }

    std::cout << "Output xmlPath:     " << xmlPath << std::endl;
    std::cout << "Output histPath:    " << histPath << std::endl;
    std::cout << "Output file prefix: " << outputPrefix << std::endl;

    if (myFile.contains("parameters","nEvents")) {
        numEvents = myFile.getInt("parameters", "numEvents");
    }

    std::cout << "Maximum number of events to process: " << numEvents << std::endl;

    BadStripsCalib* r = new BadStripsCalib(
        myFile,
        digiChain, 0, 0, 
        const_cast<char*>(outputPrefix.c_str()),
        const_cast<char*>(xmlPath.c_str()), 
        const_cast<char*>(histPath.c_str())
    );  
    r->Go(numEvents);
    r->Finish();
    r->WriteHist();
    delete r;
   
    if (attended) {
        char istop;
        std::cout << "Hit Enter key to exit: ";
        std::cin.get(istop);
    }
    return 0;
}
