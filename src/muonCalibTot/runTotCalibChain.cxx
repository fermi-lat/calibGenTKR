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

  if(argn > 1) {
    inputFile.open(argc[1]);
  }
  else {
    inputFile.open("../src/muonCalibTot/totCalibChain_option.dat");
  }

  // first line is a list of recon root files, separated by " " 
  std::string line;
  do{ getline(inputFile, line);
  } while( line[0] == '#' ); // skip the line with #.

  std::vector<std::string> digiFileNames;

  parseFileNames(digiFileNames, line);

  TChain* digiChain = new TChain("Digi");

  for(std::vector<std::string>::const_iterator itr = digiFileNames.begin();
      itr != digiFileNames.end(); ++itr) {

    std::cout << "digi file: " << *itr << endl;

    digiChain->Add(itr->c_str());

  }

  // second line is a list of recon root files, separated by " "
  do{ getline(inputFile, line);
  } while( line[0] == '#' ); // skip the line with #.

  std::vector<std::string> reconFileNames;

  parseFileNames(reconFileNames, line);

  TChain* reconChain = new TChain("Recon");

  for(std::vector<std::string>::const_iterator itr = reconFileNames.begin();
      itr != reconFileNames.end(); ++itr) {

    std::cout << "recon file: " << *itr << endl;

    reconChain->Add(itr->c_str());

  }

  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string txtFileName = line;

  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string xmlFileName = line;

  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string rootFileName = line;

  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string totConvDir = line;

  do{ getline(inputFile, line);
  } while( line[0] == '#' );
  std::string totConvRunId = line;

  totCalib calib;
  if( !calib.readTotConvFile( totConvDir.c_str(), totConvRunId.c_str() ) )
    return 1;

  if( !calib.setOutputFiles( txtFileName.c_str(), xmlFileName.c_str(), 
			     rootFileName.c_str() ) ) return 1;

  int nEvents = calib.setInputRootFiles( digiChain, reconChain );

  calib.calibChargeScale( nEvents );

}
