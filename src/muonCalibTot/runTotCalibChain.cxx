#include <fstream>
#include <iostream>
#include <string>
#include "totCalib.h"

using std::string;
using std::cout;
using std::endl;

void parseFileNames(std::vector<std::string>& fileNames, 
		    const std::string& line) 
{
  std::string::size_type pos = 0;
  for( ; ; ) {

    string::size_type i = line.find(' ', pos);
    if(i != string::npos) {
      fileNames.push_back(line.substr(pos, i-pos));
    }
    else {

      std::string lastFile = line.substr(pos);

      //make sure it is a root file name
      if(lastFile.find("root") != string::npos) {
	fileNames.push_back(lastFile);
	break;
      }
      else {
	break;
      }

    }

    pos = i + 1;
  }
}

int main(int argn, char** argc) {

  std::ifstream inputFile;
  int maxEvents = -1;

  if(argn > 1) {
    maxEvents = atoi( argc[1] );
  }
  if(argn > 2) {
    inputFile.open(argc[2]);
  }
  else {
    inputFile.open("../src/muonCalibTot/totCalibChain_option.dat");
  }

  std::string line;

  // top directory for input root files
  do{ getline(inputFile, line);
  } while( line[0] == '#' ); // skip the line with #.
  std::string rootDir = line;

  // directory name for recon root files
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string reconDir = line;

  // top directory for input root files
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string reportDir = line;

  // run ids for input root files
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::vector<std::string> runIds;
  parseFileNames( runIds, line );

  // top directory for TOT conversion factors
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string totConvDir = line;

  // run Id for TOT conversion factors
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string totConvRunId = line;

  // output directory name
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string outputDir = line;

  //dtd name for output xml
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string dtd = line;

  // muon MIP calibration or bad strip
  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string analysisType = line;

  totCalib calib( analysisType );

  if( !calib.setOutputFiles( outputDir.c_str() ) ) return 1;

  if( !calib.readTotConvXmlFile( totConvDir.c_str(), totConvRunId.c_str() ) )
    if( !calib.readTotConvFile( totConvDir.c_str(), totConvRunId.c_str() ) )
      return 1;

  if( dtd.size() > 5 ) calib.setDtd( dtd );

  int nEvents = calib.setInputRootFiles( rootDir.c_str(), reconDir.c_str(), 
					 runIds );
  if( maxEvents > 0 ) nEvents = maxEvents;

  if( !calib.readRcReports( reportDir.c_str(), runIds ) ) return 1;

  calib.calibChargeScale( nEvents );

}
