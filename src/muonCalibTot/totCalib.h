#ifndef totCalib_Class
#define totCalib_Class

#include "xml/XmlParser.h" 
#include <map>
#include <string>
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
#include "xml/Dom.h"
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

#include "xml/Dom.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "facilities/Util.h"

#include <string>
#include <iostream>
#include <fstream>

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

class totCalib {
 public:

  totCalib();
  ~totCalib();

  int setInputRootFiles( const char* digi, const char* recon );
  int setInputRootFiles( const char* rootDir, const char* reconDir,  
			 const std::vector<std::string>& runIds );
  int setInputRootFiles( TChain* digi, TChain* recon );

  bool setOutputFiles( const char* txt, const char* xml, const char* root);

  bool readTotConvFile(const char* dir, const char* runid);
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

  void retrieveCluster();

  void fillXml();//takuya

  int findTot(int planeId, TkrCluster::view viewId, int stripId);

  bool readTotConv(int layer, int view, const char* file);
  bool readTotConv(const char* file);

  bool getParam(const DOMElement* totElement,int layer,int view);

  float calcCharge(int planeId, TkrCluster::view viewId, int iStrip, 
		   int tot) const;

  static const int g_nLayer = 18;
  // g_nPlane=0 refers to top biLayer while g_nLayer=0 refers to bottom biLayer
  static const int g_nPlane = 18;
  static const int g_nFecd = 24;
  static const int g_nView = 2;
  static const int g_nStrip = 1536;

  // group of strips so that each group has enough statistics
  static const int g_nDiv = 24; //used to be 64;takuya

  TGraphErrors* m_totStrip[g_nPlane][g_nView];
  TGraphErrors* m_chargeStrip[g_nPlane][g_nView];

  TH1F* m_totHist[g_nPlane][g_nView][g_nDiv];
  TH1F* m_chargeHist[g_nPlane][g_nView][g_nDiv];

  float m_chargeScale[g_nPlane][g_nView][g_nDiv];

  TNtuple* m_tuple;

  TFile* m_reconFile;
  TTree* m_reconTree;
  ReconEvent* m_reconEvent;

  TFile* m_digiFile;
  TTree* m_digiTree;
  DigiEvent* m_digiEvent;

  TFile* m_totFile;

  // reconstructed event vertex and direction
  TVector3 m_pos, m_dir;

  int m_totX[g_nPlane][2];
  int m_totY[g_nPlane][2];

  // Index of cluster hit in the TkrSiClusters
  std::map<int, TkrCluster*> m_cluster;
 
  //file stream for log output
  std::ofstream m_log;
  // file name for output xml file
  std::string m_xmlOutput;//takuya

  float m_totQuadra[g_nLayer][g_nView][g_nStrip];
  float m_totGain[g_nLayer][g_nView][g_nStrip];
  float m_totOffset[g_nLayer][g_nView][g_nStrip];

  //xml related parameters
  int m_tower_row, m_tower_col, m_first_run, m_last_run;
  std::string m_tower_serial, m_version, m_tot_runid, m_timeStamp, 
    m_startTime, m_stopTime;

};
Double_t langaufun(Double_t *x, Double_t *par);//takuya0122
#endif
