#ifndef totCalib_Class
#define totCalib_Class

#include <map>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TStyle.h"

#include "digiRootData/DigiEvent.h"
#include "reconRootData/ReconEvent.h"
#include "IndexedVector.h"

namespace {
    int g_nTower = 18;
    int g_nLayer = 18;
    int g_nPlane = 36;
    //const int g_nFecd = 24;
    const int g_nView = 2;
    const int g_nStrip = 1536;

    // group of strips so that each group has enough statistics
    const int g_nDiv = 24; // corresponds to one front-end chip
}

class totCalib {
 public:

  totCalib();
  ~totCalib();

  void genTot(const char* digi, const char* recon, 
             const char* outputTxt, const char* outputRoot);
  void genTot(TChain* digi, TChain* recon, 
             const char* outputTxt, const char* outputRoot);

 private:

  void analyzeEvent(int nEvent);

  bool passCut(TkrRecon* tkrRecon);

  void fillTot(TkrRecon* tkrRecon); 

  void fitTot();

  TH1F* m_peaks;
  TH1F* m_widths;
  TH1F* m_chisqs;

  IndexedVector <TGraphErrors*> m_totStrip;

  IndexedVector<TH1F*> m_totHist;

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
 
  //name of ascii output tot file
  std::string m_txtOutput;
};

#endif
