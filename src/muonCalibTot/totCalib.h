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

static const float sumThrPerEvent = 60.0 / 2.0E6;
static const float occThrPerEvent = 10.0 / 1.6E6;
static const float poissonThreshold = -5.0;



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
    void setDtd( const std::string &dtd ){ m_dtd = dtd;
    std::cout << "DTD file: " << m_dtd << std::endl;
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

    void retrieveCluster();

    void fillXml();//takuya

    int findTot(int tower,int planeId, int view, int stripId);

    bool getParam(const DOMElement* totElement, int tower, int layer, int view);
    bool checkFile( const std::string );

    float calcCharge(int tower, int layer, int view, int iStrip, int tot) const;

    static const int g_nLayer = 18;
    // g_nPlane=0 refers to top biLayer while g_nLayer=0 refers to bottom biLayer
    static const int g_nPlane = 18;
    static const int g_nFecd = 24;
    static const int g_nView = 2;
    static const int g_nStrip = 1536;

    // group of strips so that each group has enough statistics
    static const int g_nDiv = 24; //used to be 64;takuya

    static const int g_nTower = 2;

    TGraphErrors* m_totStrip[g_nTower][g_nLayer][g_nView];
    TGraphErrors* m_chargeStrip[g_nTower][g_nLayer][g_nView];

    TH1F* m_totHist[g_nTower][g_nLayer][g_nView][g_nDiv];
    TH1F* m_chargeHist[g_nTower][g_nLayer][g_nView][g_nDiv];

    float m_chargeScale[g_nTower][g_nLayer][g_nView][g_nDiv];

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

    int m_tot[g_nTower][g_nLayer][g_nView][2];
    int m_lastRC0Strip[g_nTower][g_nLayer][g_nView];

    // Index of cluster hit in the TkrSiClusters
    std::map<int, TkrCluster*> m_cluster;

    //file stream for log output
    std::ofstream m_log;
    // output directory
    std::string m_outputDir;
    // dtd file
    std::string m_dtd;

    // tower list
    std::vector<int> m_towerList;

    float m_totQuadra[g_nTower][g_nLayer][g_nView][g_nStrip];
    float m_totGain[g_nTower][g_nLayer][g_nView][g_nStrip];
    float m_totOffset[g_nTower][g_nLayer][g_nView][g_nStrip];

    //xml related parameters
    int m_first_run, m_last_run;
    std::string m_tower_serial[g_nTower], m_tot_runid[g_nTower], 
      m_version, m_tag, m_dateStamp, m_timeStamp, m_startTime, m_stopTime;

    // bad strips analysis related stuff
    static const int g_nWafer = 4, g_nBad=6;
    void fillOccupancy();
    void findBadStrips( int );
    void openBadStripsXml( std::ofstream &xmlFile, std::ifstream &dtd ); 
    void fillBadStrips();
    void fillTowerBadStrips( std::ofstream &xmlFile, const int tower, 
			     const int nBad=g_nBad-1 );

    bool m_badStrips;
    std::vector<int> m_deadStrips[g_nTower][g_nLayer][g_nView][g_nBad];
    TH1F* m_nHits[g_nTower][g_nLayer][g_nView][g_nWafer];
    TH1F* m_aPos[g_nWafer];
    TH1F* m_occDist, *m_poissonDist;

};

Double_t langaufun(Double_t *x, Double_t *par);//takuya0122
#endif
