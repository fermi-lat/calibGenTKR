#include <fstream>
#include <iostream>
#include <string>
#include "totCalib.h"

#include "facilities/commonUtilities.h"

int main(int argn, char** argc) {
  facilities::commonUtilities::setupEnvironment();

  int maxEvents = -1;
  std::string jobXml = "../src/muonCalibTot/calibTkrJobOptions.xml";
  std::string defJob = "None";
  
  std::cout << "Start doMuonCalibTot" << std::endl;
  
  if(argn > 1) {
    maxEvents = atoi( argc[1] );
    std::cout << "# of events spesified: " << maxEvents << std::endl;
  }
  if(argn > 2) {
    defJob = argc[2];
    std::cout << "job option spesified: " << defJob << std::endl;
  }
  if(argn > 3) {
    jobXml = argc[3];
    std::cout << "job option file spesified: " << jobXml << std::endl;
  }
  
  totCalib calib( jobXml, defJob );
  calib.analyze( maxEvents );
  
}

