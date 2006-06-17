#ifndef totCalib_Class
#define totCalib_Class

#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMCharacterData.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMTreeWalker.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"

#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "facilities/Util.h"

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "calibTkrUtil/TkrHits.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

static const float poissonThreshold = -5.0;
static const float minEff = 0.85;
const float maxProbT = -5.0;
static const float maxProb = -3.5, maxProbW = -3.5, maxContSum = -8.0, maxProbSum = -10.0;
static const float probThreshold = -5.0;

const float maxChisq = 1.75;
const float maxFracErr = 0.015;
const float minFracErr = 0.005;

const float maxTot = 20.0;

//
// class poisson
//
class poissonFunc{
 public:
  poissonFunc();
  ~poissonFunc(){;};

  double getProb( double, int );
  float getLogProb( float, int );
  float getLogIntProb( double, int );
  
 private:
  double factorial[200];
  float logFactorial[200];
};

//
// ***** class totCalib *****
//
class totCalib:public TkrHits {
 public:
  
  totCalib( const std::string, const std::string );
  ~totCalib();
  
  void getTimeStamp();
  void initHists();
  void initTotHists();
  
  bool readJobOptions( const std::string, const std::string );
  void parseRunIds(std::vector<std::string>&, const std::string& );
  void splitWords(std::vector<std::string>&, const std::string& );

  int setInputRootFiles( const char* digi, const char* recon );
  int setInputRootFiles( const char* rootDir, const char* digiPrefix,
			 const char* reconPrefix,  
			 const std::vector<std::string>& runIds );
  int setInputRootFiles( TChain* digi, TChain* recon );
  bool addToChain( const char* rootDir, const char* digiPrefix,
		   const char* reconPrefix, 
		   const std::vector<std::string> &runIds,
		   TChain* digiChain, TChain* reconChain );
  
  bool setOutputFiles( const char* outputDir );
  void setDtdDir( const std::string &dtdDir ){ m_dtdDir = dtdDir;
  std::cout << "DTD directory: " << m_dtdDir << std::endl;
  };
  
  bool readInputXmlFiles( const std::string, 
			  const std::vector<std::string>& runIds );  
  bool readBadStripsXmlFile( const std::string dir, const std::string runid );
  bool readBadStripsXmlFile( const std::string filename, bool hotStrips );
  bool readTotConvXmlFile( const std::string dir, const std::string runid );
  bool readTotConvXmlFile( const std::string filename );
  bool readLatcTfeXmlFile( const std::string filename );

  bool readBadStripsTxtFiles( const std::string dir, const std::string names);
  bool readBadStripsTxtFile( const std::string filename );

  bool readInputHistFiles( const std::string, 
			   const std::vector<std::string>& );
  bool readInputHistFiles( const std::string dir, 
			   const std::string prefix, 
			   const std::vector<std::string> &runIds );
  bool readHists( TFile*, UInt_t, UInt_t );
  inline void histAdd( TH1F* hist, TFile* hfile, const char* name ){
    TH1F* h1f = (TH1F*)hfile->FindObjectAny( name );
    if( h1f ) hist->Add( h1f );
  };

  bool readRcReports( const char* reportDir, 
		      const std::vector<std::string>& runIds );
  bool readHotStrips( const char* reportDir, 
		      const std::vector<std::string>& runIds );
  
  bool parseRcReport( const char* reportFile );
  void totCalib::getDate( const char* str, std::string& date );
  
  void analyze( int );
  
 private:
  
  void analyzeEvents();
  
  void fillTot();   
  void getTot();
  void fitTot();
  
  void fillXml();//takuya
  void fillTowerChargeScales( std::ofstream &xmlFile, const int tower );
  void openChargeScaleXml( std::ofstream &xmlFile, std::ifstream &dtd, const std::string tot_runid );
  
  int findTot(int tower,int planeId, int view, int stripId);
  
  bool getParam(const DOMElement* totElement, layerId lid, std::vector<std::string> keywords );
  bool checkFile( const std::string );
  
  float calcCharge( layerId lid, int iStrip, int tot) const;
  
  std::vector<TH1F*> m_chargeHist;

  TFile* m_reconFile;
  TTree* m_reconTree;
  
  TFile* m_digiFile;
  TTree* m_digiTree;
 
  int m_tot[g_nTower][g_nLayer][g_nView][2];
  
  // output directory
  std::string m_outputDir;
  // dtd file
  std::string m_dtdDir;

  std::string m_nameType;
  
  //xml related parameters
  int m_first_run, m_last_run;
  std::string m_startTime, m_stopTime, m_dateStamp, m_timeStamp;;
  
  // bad strips analysis related stuff
  void findBadStrips();
  void openBadStripsXml( std::ofstream &xmlFile, std::ifstream &dtd ); 
  void fillBadStrips();
  void fillTowerBadStrips( std::ofstream &xmlFile, const int tower, 
			   const bool saveHotStrips=false,
			   const int nBad=g_nBad-1 );
  void combineBadChannels( layerId );
  void fixedDisp( int, bool=true );
  void calculateEfficiency();
  
  bool m_correctedTot, m_histMode;
  float m_totAngleCF, m_totThreshold, m_totGain, m_totQuad, m_RSigma, m_GFrac, m_peakMIP;
  
};
// Gaussian convolved Landau function
Double_t langaufun(Double_t *x, Double_t *par); 
// Two Gaussian convolved Landau function
Double_t langau2fun(Double_t *x, Double_t *par); 
#endif
