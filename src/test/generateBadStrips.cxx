#include "badStripsCalib.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>
#include "xml/IFile.h"

/* to do:
  figure out better local average, for now we use simple average
  pass string for output files
  input params from text file?
  */

/** @file generateBadStrips
@brief Driver for bad strips calibration
*/

std::string stripBlanks(std::string &str) {
    std::string temp = str;
    if (temp.size()==0) return temp;
    std::string::size_type pos;
    pos = temp.find_first_not_of(" ", 0);
    temp = temp.substr(pos);
    pos = temp.find_last_not_of(" ", string::npos);
    temp = temp.substr(0,pos+1);
    return temp;
}

int splitString(std::string &input, std::string &LH, std::string &RH, char* delim) {
    std::string::size_type pos;
    pos = input.find(delim);
    if (pos!=string::npos) {
        LH = input.substr(0,pos);
        RH = input.substr(pos);
    } else {
        LH = input;
        RH = "";
    }
    LH = stripBlanks(LH);
    RH = stripBlanks(RH);
    return pos;
}

int main(int argn, char** argc) {
    
#ifdef WIN32
    gSystem->Load("libTree.dll");   
#endif

    bool attended = false;
    std::string temp;
 
    std::string sourceFilePath;
    std::string sourceFile;
    std::string outputString;
    std::string path = ::getenv("CALIBGENTKRROOT");
    unsigned int numEvents = 5000000;
 
    std::string xmlFile(path+"/src/test/options.xml");

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
        temp = myFile.getString("parameters", "sourceFileList");
        sourceFileString = stripBlanks(temp);
    }

    int nFiles = 0;
    std::cout << "Input files:" << std::endl;
    TChain* digiChain = new TChain("Digi");    
    std::string::size_type pos;
    temp = sourceFileString;
    while(temp!="") {
        pos = splitString(temp, sourceFile, temp, " ");
        if (sourceFile!="") nFiles++;
        digiChain->Add((sourceFilePath+sourceFile).c_str());
        std::cout << "   " << nFiles << ") " << sourceFile << std::endl;
    }

    if (myFile.contains("parameters","outputString")) {
        temp = myFile.getString("parameters", "outputString");
        outputString = temp;
    }

    std::cout << "Output file prefix in directory /output: " << outputString << std::endl;

    if (myFile.contains("parameters","numEvents")) {
        numEvents = myFile.getInt("parameters", "numEvents");
    }

    std::cout << "Maximum number of events to process: " << numEvents << std::endl;

    std::string outputPrefix(path+"/output/"+outputString);

    BadStripsCalib* r = 
        new BadStripsCalib(digiChain, 0, 0, const_cast<char*>(outputPrefix.c_str()));  
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







