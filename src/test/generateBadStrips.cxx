#include "badStripsCalib.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <string>

/** @file generateBadStrips.cxx
@brief Driver for bad strips calibration
*/


int main(int argn, char** argc) {
    
#ifdef WIN32
    gSystem->Load("libTree.dll");   
#endif

    bool attended = false;
    string temp;
 
    string sourceFilePath;
    string sourceFile;
    string path = ::getenv("CALIBGENTKRROOT");
    unsigned int numEvents = 5000000;

    string xmlPath(path+"/output/");
    string histPath(path+"/output/");
    string outputPrefix("test");
 
    string xmlFile(path+"/src/test/options.xml");

    if(argn > 1) {
        xmlFile = argc[1];
        cout << "Reading in user-specified options file: " <<  xmlFile 
            << endl << endl;
    }
    
    xml::IFile myFile(xmlFile.c_str());

    if (myFile.contains("parameters","sourceFilePath")) {
        temp = myFile.getString("parameters", "sourceFilePath");
        sourceFilePath = stripBlanks(temp);
    }

    cout << "Sourcefile path: " << sourceFilePath << endl;

    string sourceFileString;
    if (myFile.contains("parameters","sourceFileList")) {
        temp = myFile.getString("parameters", "sourceFileList");
        sourceFileString = stripBlanks(temp);
    }

    int nFiles = 0;
    cout << "Input files:" << endl;
    TChain* digiChain = new TChain("Digi");    
    string::size_type pos;
    temp = sourceFileString;
    while(temp!="") {
        pos = splitString(temp, sourceFile, temp, " ");
        if (sourceFile!="") nFiles++;
        digiChain->Add((sourceFilePath+sourceFile).c_str());
        cout << "   " << nFiles << ") " << sourceFile << endl;
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

    cout << "Output xmlPath:     " << xmlPath << endl;
    cout << "Output histPath:    " << histPath << endl;
    cout << "Output file prefix: " << outputPrefix << endl;

    if (myFile.contains("parameters","numEvents")) {
        numEvents = myFile.getInt("parameters", "numEvents");
    }

    cout << "Maximum number of events to process: " << numEvents << endl;

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
        cout << "Hit Enter key to exit: ";
        cin.get(istop);
    }
    return 0;
}







