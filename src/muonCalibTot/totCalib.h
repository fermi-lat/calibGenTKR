#ifndef totCalib_Class
#define totCalib_Class

#ifdef OLD_RECON
#include "xml/Dom.h"
#else
#include "xmlBase/Dom.h"
#endif
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
#ifdef OLD_RECON
#include "xml/XmlParser.h"
#include "xml/Dom.h"
#else
#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#endif

#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "facilities/Util.h"

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "digiRootData/DigiEvent.h"
#include "reconRootData/ReconEvent.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

//static const float sumThrPerEvent = 60.0 / 2.0E6;
//static const float occThrPerEvent = 10.0 / 1.6E6;
static const float poissonThreshold = -5.0;
static const float minEff = 0.85;
const float minProb = -5.0;
static const float maxProb = -2.5, maxProbW = -3.0, maxProbT = -5.0, maxProbSum = -7.0;
static const float probThreshold = -5.0;

const float maxChisq = 1.6;
const float maxFracErr = 0.015;

const int nTotHistBin = 200;
const float maxTot = 20.0;

static const float stripPitch = 0.228; // strip pitch

static const int g_nLayer = 18;
// g_nPlane=0 refers to top biLayer while g_nLayer=0 refers to bottom biLayer
#ifdef OLD_RECON
static const int g_nPlane = 18;
#endif
static const int g_nView = 2;
static const int g_nUniPlane = g_nLayer*g_nView;

static const int g_nFecd = 24;
static const int g_nStrip = 1536;
static const int g_nTower = 8;

// group of strips so that each group has enough statistics
static const int g_nDiv = 24; //used to be 64;takuya
  
static const int g_nWafer = 4, g_nBad=6, g_nTime=5, g_nMerge=4;

const float ladderGap = 2.148;
const float posZ[g_nView][g_nLayer] = { 
  {44.765, 74.335, 108.965, 139.538, 175.165, 205.738, 241.365,  271.14, 305.965, 335.74, 370.565, 400.34, 435.165, 464.94, 499.765, 529.54, 564.365, 594.14},
  {42.315, 76.865, 106.435, 142.065, 172.638, 208.265, 238.838, 273.665, 303.44, 338.265, 368.04, 402.865, 432.64, 467.465, 497.24, 532.065, 561.84, 596.665} };

//
// ****** class layerId *****
//
class layerId{
 public:
      layerId(){;};
      layerId( int, int, int twr=0 );
      layerId( int, std::string, int twr=0 );
      layerId( int );
      ~layerId(){;};
      
      void setLayer( int, int, int twr=0 );
      void setTray( int, std::string, int twr=0 );
      void setUniPlane( int, int twr=0 );
      void setTower( int tw ){ tower = tw; };
      
      void trayToUniPlane();
      void trayToLayer();
      void layerToTray();
      inline void layerToUniPlane(){ layerToTray(); trayToUniPlane(); };
      void uniPlaneToTray();
      inline void uniPlaneToLayer(){ uniPlaneToTray(); trayToLayer(); };
      
      int tower;
      int layer;
      int view;
      int uniPlane;
      int tray;
      std::string which;
      
};

//
// cluster class
//
class Cluster{
 public:
  Cluster( int, int );
  ~Cluster(){;};
  bool addStrip( int );

 private:
  int tower;
  int uniPlane;
  int firstStrip; // first strip number
  int lastStrip; // last strip number
  int tot;
  TVector3 pos; // position

 public:
  inline void setId( int tw, int unp ){ tower = tw; uniPlane = unp; };
  inline int getTowerId(){ return tower; };
  inline int getUniPlane(){ return uniPlane; };
  inline int getFirstStrip(){ return firstStrip; };
  inline int getLastStrip(){ return lastStrip; };
  inline int getSize(){ return lastStrip-firstStrip+1; };
  inline void setTot( int val ){ tot = val; };
  inline int getRawToT(){ return tot; };
  inline TVector3 getPosition(){ return pos; };
  inline void setXYZ( float x, float y, float z ){ pos.SetXYZ( x, y, z ); };
};

//
// bad strips variables
//
struct badStripVar{
  int tHits[g_nStrip];
  int lHits[g_nStrip];
  int nHits[g_nStrip][g_nWafer][g_nTime]; 
  std::vector<int> badStrips[g_nBad];
};

//
// tot calibration variable
//
struct totCalibVar{
  float totQuadra[g_nStrip];
  float totGain[g_nStrip];
  float totOffset[g_nStrip];
  float chargeScale[g_nDiv];
};

//
//
// tower variable class
//
class towerVar{
 public:
  towerVar( int, bool );
  ~towerVar(){;};
  void saveHists();

  int towerId;
  float center[g_nView];
  std::string hwserial, runid;

  // clusters
  std::vector<Cluster> digiClusters[g_nUniPlane];
  const TkrCluster* reconClusters[g_nUniPlane];

  // do not use histogram since it will slow down significantly.
  // these will be converted to histograms later in saveHits method.
  int rHits[g_nUniPlane][g_nStrip];
  int dHits[g_nUniPlane][g_nStrip];
  std::vector<badStripVar> bsVar;
  std::vector<totCalibVar> tcVar;
};

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
class totCalib {
 public:
  
  totCalib( const std::string );
  ~totCalib();
  
  int setInputRootFiles( const char* digi, const char* recon );
  int setInputRootFiles( const char* rootDir, const char* digiPrefix,
			 const char* reconPrefix,  
			 const std::vector<std::string>& runIds );
  int setInputRootFiles( TChain* digi, TChain* recon );
  
  bool setOutputFiles( const char* outputDir );
  void setDtdDir( const std::string &dtdDir ){ m_dtdDir = dtdDir;
  std::cout << "DTD directory: " << m_dtdDir << std::endl;
  };
  
  bool readTotConvXmlFile(const char* dir, const char* runid);
  bool readRcReports( const char* reportDir, 
		      const std::vector<std::string>& runIds );
  
  bool parseRcReport( const char* reportFile );
  void totCalib::getDate( const char* str, std::string& date );
  
  void calibChargeScale( int );
  
 private:
  
  void analyzeEvent(int nEvent);
  
  bool passCut();
  
  void fillTot(); 
  
  void getTot(); 
  
  void fitTot();
  
  void getDigiClusters();
  void getReconClusters();
  void selectGoodClusters();
  layerId getLayerId( const TkrCluster* );
  layerId getLayerId( Cluster* );
  bool closeToTrack( const TkrCluster*, TkrCluster*[g_nTower][g_nUniPlane] );
  
  void fillXml();//takuya
  void fillTowerChargeScales( std::ofstream &xmlFile, const int tower );
  void openChargeScaleXml( std::ofstream &xmlFile, std::ifstream &dtd, const std::string tot_runid );
  
  int findTot(int tower,int planeId, int view, int stripId);
  
  bool getParam(const DOMElement* totElement, layerId lid );
  bool checkFile( const std::string );
  
  float calcCharge( layerId lid, int iStrip, int tot) const;
  
  TH1F *m_nTrackDist, *m_maxHitDist, *m_numClsDist, *m_dirzDist, 
    *m_armsDist, *m_brmsDist[g_nLayer/3];
  
#ifdef FULLHIST
  TGraphErrors* m_totStrip[g_nTower][g_nLayer][g_nView];
  TGraphErrors* m_chargeStrip[g_nTower][g_nLayer][g_nView];
  
  TH1F* m_totHist[g_nTower][g_nLayer][g_nView][g_nDiv];
#endif
  TH1F* m_chargeHist[g_nTower][g_nUniPlane][g_nDiv];
  TH1F *m_fracErrDist, *m_chisqDist, *m_fracBatTot;
  //float m_chargeScale[g_nTower][g_nLayer][g_nView][g_nDiv];

  TNtuple* m_tuple;
  
  TFile* m_reconFile;
  TTree* m_reconTree;
  ReconEvent* m_reconEvent;
  
  TFile* m_digiFile;
  TTree* m_digiTree;
  DigiEvent* m_digiEvent;
  
  TFile* m_totFile; 
 
  // reconstructed track and direction
#ifdef OLD_RECON
  TkrKalFitTrack* m_track;
#else
  TkrTrack* m_track;
#endif
  TVector3 m_pos, m_dir;
  
  std::vector<Cluster*> m_clusters;

  int m_tot[g_nTower][g_nLayer][g_nView][2];
  int m_lastRC0Strip[g_nTower][g_nLayer][g_nView];
  
  //file stream for log output
  std::ofstream m_log;
  // output directory
  std::string m_outputDir;
  // dtd file
  std::string m_dtdDir;
  
  // tower variables
  std::vector<int> m_towerList;
  int m_towerPtr[g_nTower];
  std::vector<towerVar> m_towerVar;
  std::vector<int> m_trackTowerList;
  
  //float m_totQuadra[g_nTower][g_nLayer][g_nView][g_nStrip];
  //float m_totGain[g_nTower][g_nLayer][g_nView][g_nStrip];
  //float m_totOffset[g_nTower][g_nLayer][g_nView][g_nStrip];
  
  //xml related parameters
  int m_first_run, m_last_run;
  //std::string m_tower_serial[g_nTower], m_tot_runid[g_nTower];
  std::string m_version, m_tag, m_dateStamp, m_timeStamp, m_startTime, m_stopTime;
  
  // bad strips analysis related stuff
  void fillOccupancy( int );
  void findBadStrips();
  void openBadStripsXml( std::ofstream &xmlFile, std::ifstream &dtd ); 
  void fillBadStrips();
  void fillTowerBadStrips( std::ofstream &xmlFile, const int tower, 
			   const int nBad=g_nBad-1 );
  void fixedDisp( int, bool=true );
  
  bool m_badStrips;
  TH1F* m_aPos[g_nWafer+1];
  TH1F *m_occDist, *m_poissonDist, *m_lrec, *m_ldigi, *m_lcls, *m_locc, *m_leff, *m_ltrk, *m_dist;
  
};

Double_t langaufun(Double_t *x, Double_t *par);//takuya0122
#endif
