/** @file badStripsCalib.h
@brief header file for bad strips calibration
*/
/** @class BadStripsCalib
* @brief   This class constains the code to do the bad strips calibration.
*
* It is derived from the example RootTreeAnalysis
*/


#ifndef BadStripsCalib_h
#define BadStripsCalib_h 1

#if !defined(__CINT__)
// Need these includes if we wish to compile this code
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCollection.h"  // Declares TIter
#include "TStopwatch.h"
#include "digiRootData/DigiEvent.h"
//#include "reconRootData/ReconEvent.h"
//#include "mcRootData/McEvent.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#else  // for interactive use
#include "iostream.h"
#include "string.h"
class DigiEvent;
class ReconEvent;
class McEvent;
#endif

#include "xml/IFile.h"

namespace {
    string stripBlanks(string &str) {
        // strips off leading and trailing blanks
        string temp = str;
        if (temp.size()==0) return temp;
        string::size_type pos;
        pos = temp.find_first_not_of(" ", 0);
        temp = temp.substr(pos);
        pos = temp.find_last_not_of(" ", string::npos);
        temp = temp.substr(0,pos+1);
        return temp;
    }

    int splitString(string &input, string &LH, string &RH, char* delim) {
        // splits off leftmost token from delim-delimited string
        string::size_type pos;
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
}


class BadStripsCalib {
public :
    /// specific to this application
    /// number of planes per tower for this detector
    int m_nPlanes;
    /// list of tower numbers (also gives number of towers)
    std::vector<int> m_towerNums;
    /// hist Ids
    int m_histId[16];
    /// maximum occupancy parameter
    double m_maxOccupancy;
    /// vector of histograms
    std::vector<TH1F*> m_tkrHists;

    /// Histogram file
    TFile       *histFile;
    /// Input digitization file
    TFile       *digiFile;   
    /// Input reconstruction file
    TFile       *reconFile;  
    /// Input monte carlo file
    TFile       *mcFile;     
    /// pointer to the digi tree
    TTree       *digiTree;    
    /// pointer to the reconstruction tree
    TTree       *reconTree;
    /// pointer to the monte carlo tree
    TTree       *mcTree;      
    /// Optional TChain input
    TChain      *m_digiChain, *m_recChain, *m_mcChain;
    /// pointer to a DigiEvent
    DigiEvent   *evt;
    /// pointer to a ReconEvent
    //ReconEvent  *rec;
    /// Pointer to a McEvent
    //McEvent     *mc;
    /// name of the output histogram ROOT file
    char        *m_histFileName;
    /// xml output path
    char        *m_xmlPath;
    /// hist output path
    char        *m_histPath;

    /// output prefix
    std::string  m_prefix;
    /// Arrays that contain pointers to the TFile, TTree, and TChains
    TObjArray   *fileArr, *treeArr, *chainArr;
    
	/// Default ctor, requires that that user calls BadStripsCalib::Init
	/// to setup access to specific ROOT files.
    BadStripsCalib(); 
    
	/// Standard ctor, where user provides the names of the input root files
	/// and optionally the name of the output ROOT histogram file
    BadStripsCalib(
        xml::IFile& myFile,
        const char* digiFileName, 
        const char* reconFileName="", 
        const char* mcFileName="",
        char* prefix="test",
        char* xmlPath="",
        char* histPath=""
    );
       

	/// Special ctor which accepts TChains for input files
    BadStripsCalib( 
        xml::IFile& myFile,
        TChain *digiChain, 
        TChain *recChain = 0, 
        TChain *mcChain = 0, 
        char* prefix = "",
        char* xmlPath = "",
        char* histPath = ""
    );

    ~BadStripsCalib();  

    /// start next Go with this event
    void StartWithEvent(Int_t event) { m_StartEvent = event; };  
    /// reset for next Go to start at beginning of file 
    void Rewind() { m_StartEvent = 0; }; 
    /// allow the user to specify their own file name for the output ROOT file
    void SetHistFileName(char *histFileName) { 
        m_histFileName = histFileName;
    }
    /// re-init with these input Root files
    void Init(  const char* digiFileName="", 
        const char* reconFileName="", 
        const char* mcFileName=""); 
    /// define user histograms, ntuples and other output objects that will be saved to output
    void HistDefine();   
    /// make list of user histograms and all objects created for output
    void MakeHistList(); 
    /// write the existing histograms and ntuples out to file
    void WriteHist() { if (histFile) histFile->Write(); }; 
    /// Reset() all user histograms
    void HistClear(); 
    /// Retrieve a pointer to an object stored in our output ROOT file
    TObject* GetObjectPtr(const char *tag) { return (m_histList->FindObject(tag)); };
    /// process events
    void Go(Int_t numEvents=100000); 
    /// returns number of events in all open files
    UInt_t GetEntries() const;
    /// retrieve a pointer to event number.
    UInt_t GetEvent(UInt_t ievt);
    /// finish up the processing
    void Finish();
    
private:
    /// starting event number
    Int_t m_StartEvent;
    /// list of user histograms
    THashList *m_histList;
        
    /// reset all member variables
    void Clear(); 

    /// read in relevant params and set up
    void GetOptions(xml::IFile& myFile);
    /// Setup the Monte Calro output histograms
    void McHistDefine();
	/// Setup the Digitization output histograms
    void DigiHistDefine();
	/// Setup the Reconstruction output histograms
    void ReconHistDefine();
    
    /// event processing for the monte carlo data
    void McData();

    /// event processing for the digi TKR data
    void DigiTkr();
	/// event processing for digi CAL data
    void DigiCal();
	/// event processing for digi ACD data
    void DigiAcd();
    
    /// event processing for the recon TKR data
    void ReconTkr();
	/// event processing for the recon CAL data
    void ReconCal();
	/// event processing for the recon ACD data
    void ReconAcd();
    
};


inline BadStripsCalib::BadStripsCalib()
{
    Clear();
}

inline BadStripsCalib::BadStripsCalib(
                                   xml::IFile& myFile,
                                   const char* digiFileName, 
                                   const char* reconFileName, 
                                   const char* mcFileName,
                                   char* prefix,
                                   char* xmlPath,
                                   char* histPath
    )
    :m_prefix(prefix), m_xmlPath(xmlPath), m_histPath(histPath)

{
	// Purpose and Method:  Standard constructor where the user provides the 
	//  names of input ROOT files and optionally the name of the output ROOT
	//  histogram file.

    Clear();

    //m_prefix   =  prefix;
    //m_xmlPath  = xmlPath;
    //m_histPath = histPath;
    printf(" opening files:\n\tdigi:\t%s\n\trecon:\t%s\n\tmc:\t%s\n",
		digiFileName, reconFileName, mcFileName);
    GetOptions(myFile);
      
    std::string histName;
    histName =  m_prefix+"_hist.root";
    std::cout << "Histogram file will be: " << histName << std::endl;
    SetHistFileName(const_cast<char*>(histName.c_str()));
    HistDefine();
    MakeHistList();
    
    Init(digiFileName, reconFileName, mcFileName);
}

inline BadStripsCalib::BadStripsCalib(
                                      xml::IFile& myFile,
                                      TChain *digiChain, 
                                      TChain *recChain, 
                                      TChain *mcChain, 
                                      char* prefix,
                                      char* xmlPath,
                                      char* histPath
    )
    :m_prefix(prefix), m_xmlPath(xmlPath), m_histPath(histPath)
{
    Clear();
    
    //m_prefix   = prefix;
    //m_xmlPath  = xmlPath;
    //m_histPath = histPath;
    GetOptions(myFile);
    std::string histName;
    histName = m_histPath+m_prefix+"_hist.root";
    SetHistFileName(const_cast<char*>(histName.c_str()));
    HistDefine();
    MakeHistList();
    
    if (chainArr) delete chainArr;
    chainArr = new TObjArray();
    
    /*
    if (mcChain != 0) {
        m_mcChain = mcChain;
        mc = 0;
        m_mcChain->SetBranchAddress("McEvent",&mc);
        chainArr->Add(m_mcChain);
    */

    if (digiChain != 0) {
        evt = 0;
        m_digiChain = digiChain;
        m_digiChain->SetBranchAddress("DigiEvent",&evt);
        chainArr->Add(m_digiChain);
    }
    
    /*
    if (recChain != 0) {
        m_recChain = recChain;
        rec = 0;
        m_recChain->SetBranchAddress("ReconEvent",&rec);
        chainArr->Add(m_recChain);
    }
    */

    m_StartEvent = 0;   
}

inline void BadStripsCalib::GetOptions(xml::IFile& myFile)
{
    int i;
    m_nPlanes = 0;
    bool overridden = false;
    std::string temp1 = "";
    if (myFile.contains("parameters","detectorType")) {
        std::string temp = myFile.getString("parameters", "detectorType");
        temp = stripBlanks(temp);
        if(temp=="") {temp = "EM1";} else {temp1 = temp;}
        if(temp=="EM1") {
            m_nPlanes = 8;
            m_towerNums.push_back(0);
        } else if(temp=="EM2") {
            m_nPlanes = 10;
            m_towerNums.push_back(0);
        }else if(temp=="LAT_2Towers") {
            m_nPlanes = 36;
            m_towerNums.push_back(8);
            m_towerNums.push_back(9);
        }else if(temp=="LAT_Full") {
            m_nPlanes = 36;
            int tower;
            for (tower=0; tower<16; ++tower) {
                m_towerNums.push_back(tower);
            }
        } else {
            std::cout << "no valid detector found" << std::endl;
        }
    }
    // direct specification overrides detectorType

    if (myFile.contains("parameters","nPlanes")) {
        int temp = myFile.getInt("parameters", "nPlanes");
        overridden = true;
        if(temp>0) {
            m_nPlanes = temp;
        }
    }
    
    if(myFile.contains("parameters","towerNumbers")) {
        m_towerNums.clear();
        m_towerNums = myFile.getIntVector("parameters", "towerNumbers");
        overridden = true;
    }  

    if (temp1=="") {
        std::cout << "No standard configuration chosen" << std::endl;
    } else {
        std::cout << "Configuration: " << temp1 << std::endl;
    }
    if(temp1>"" && overridden) {
        std::cout << "Standard parameters have been overridden" << std::endl;
    }
    std::cout << "number of planes: " << m_nPlanes << std::endl;
    std::cout << "Towers in config: " ;
    for (i=0;i<m_towerNums.size(); ++i) {
        std::cout << m_towerNums[i] << " " ;
    }
    std::cout << std::endl;

    m_maxOccupancy = 0.01;
    if (myFile.contains("parameters","maxOccupancy")) {
        m_maxOccupancy = myFile.getDouble("parameters", "maxOccupancy");
    }

    // -1 means a missing tower
    for (i=0; i<16; ++i) {
        m_histId[i] = -1;
    }
    // fill in histId sequence for existing towers
    for (i=0;i<(int)m_towerNums.size(); ++i) {
        m_histId[m_towerNums[i]] = i;
    }
}

inline BadStripsCalib::~BadStripsCalib() {
    histFile->Close();
    
    //if (m_histList) delete m_histList;
    
    if (histFile) delete histFile;
    
    if (digiFile) delete digiFile;
    //if (reconFile) delete reconFile;
    //if (mcFile) delete mcFile;
    
    if (evt) { 
		evt->Clear(); 
		delete evt;
	}
    /*
    if (rec) {
		rec->Clear();
		delete rec;
	}
    */
    /*
    if (mc) {
		mc->Clear();
		delete mc;
	*/
    
    digiTree = 0;
    //reconTree = 0;
    //mcTree = 0;
    
    if (fileArr) delete fileArr;
    if (treeArr) delete treeArr;
    if (chainArr) delete chainArr;

	Clear();
}

inline void BadStripsCalib::Init(const char* digiFileName, const char* reconFileName, const char* mcFileName)
{
    // Purpose and Method:  Re-initialize file, tree, event pointers, using the 
	//   input ROOT files.  Histograms are *not* cleared.
    
    if (fileArr) delete fileArr;
    fileArr = new TObjArray();
    
    if (treeArr) delete treeArr;
    treeArr = new TObjArray();
         
    /*
    if (mcFile) {
        delete mc; 
        mc = 0;
        mcTree = 0;
        delete mcFile; 
        mcFile = 0;
    }
    */
    
    /*
    std::string mcName = mcFileName;
    if (mcName != "") {
        mcFile = new TFile(mcFileName);
        if (mcFile->IsOpen() == kTRUE) {
            mcTree = (TTree*)gDirectory->Get("Mc");
            mc = 0;
            mcTree->SetBranchAddress("McEvent",&mc);
            fileArr->Add(mcFile);
            treeArr->Add(mcTree);
        } else {
            mcFile = 0;
            std::cout << "mc data file could not be opened!!" << std::endl;
        }
    }
    */

    if (digiFile) {
        delete evt; 
        evt = 0;
        digiTree = 0;
        delete digiFile; 
        digiFile = 0;
    }
    
    if (digiFileName != "") {
        digiFile = new TFile(digiFileName);
        if (digiFile->IsOpen() == kTRUE) {
            digiTree = (TTree*)gDirectory->Get("Digi");
            evt = 0;
            digiTree->SetBranchAddress("DigiEvent",&evt);
            fileArr->Add(digiFile);
            treeArr->Add(digiTree);
        } else {
            digiFile = 0;
            std::cout << "digi data file could not be opened!!" << std::endl;
        }
    }
    
    /*
    if (reconFile) {
        delete rec; 
        rec = 0;
        reconTree = 0;
        delete reconFile;
        reconFile = 0;
    }
    */
    
    /*
    std::string recName = reconFileName;
    if (recName != "") {
        std::cout << "reconFileName *" << reconFileName << "*" << std::endl;
        reconFile = new TFile(reconFileName);
        if (reconFile->IsOpen() == kTRUE) {
            reconTree = (TTree*)gDirectory->Get("Recon");
            rec = 0;
            reconTree->SetBranchAddress("ReconEvent",&rec);
            fileArr->Add(reconFile);
            treeArr->Add(reconTree);
        } else {
            reconFile = 0;
            std::cout << "recon data file could not be opened!!" << std::endl;
        }
    }
    */
    
    m_StartEvent = 0;
    
}


inline UInt_t BadStripsCalib::GetEvent(UInt_t ievt) {
    // Purpose and Method:  Get the event, ievt, for all trees
    //    We could be processing single files or chains, 
	//    This routine handles both casees.

    // if using regular trees - we check the array of open trees and
    // move the event pointer to the requested event
    UInt_t nb = 0;
    if (treeArr) {
        for (Int_t i = 0; i < treeArr->GetEntries(); i++) {
            nb += ((TTree*)treeArr->At(i))->GetEvent(ievt);
        }
        return nb;
    }
    
    // if using chains, check the array of chains and move
    // the event pointer to the requested event
    if (chainArr) {
        for (Int_t i = 0; i < chainArr->GetEntries(); i++) {
            nb += ((TChain*)chainArr->At(i))->GetEvent(ievt);
        }
        return nb;
    }

  return nb;
}


inline UInt_t BadStripsCalib::GetEntries() const {    
    // Purpose and Method:  Determine the number of events to iterate over
    //   checking to be sure that the requested number of events is less than
    //   the min number of events in all files

    UInt_t nentries = 0;
    if (treeArr) {
        nentries = ((TTree*)treeArr->At(0))->GetEntries();
        for (Int_t i = 1; i < treeArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, (UInt_t)((TTree*)treeArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    if (chainArr) {
        nentries = ((TChain*)chainArr->At(0))->GetEntries();
        for (Int_t i = 1; i < chainArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, (UInt_t)((TChain*)chainArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    return nentries;
}


inline void BadStripsCalib::MakeHistList() {
    // Purpose and Method:  Make a THashList of histograms
    //   This avoids the need to refresh the histogram pointers
    
    if (m_histList) delete m_histList;
    
    m_histList = new THashList(30, 5);
    
    TList* list = histFile->GetList();
    TIter iter(list);
    
    TObject* obj = 0;
    
    while (obj=iter.Next()) {
        m_histList->Add(obj);
    }
}

inline void BadStripsCalib::HistClear() {
    // Purpose and Method:  Clear histograms by iterating over the THashList
    
    if (!m_histList) return;
    
    TIter iter(m_histList);
    
    TObject* obj = 0;
    
    while ( obj=(TObject*)iter.Next() ) {
        ((TH1*)obj)->Reset();        
    }
}

inline void BadStripsCalib::Clear() {
    histFile = 0; 
    m_histList = 0;
    
    digiFile = 0; 
    //reconFile = 0;
    //mcFile = 0;
    
    digiTree = 0; 
    //reconTree = 0;
    //mcTree = 0;
    
    m_digiChain = 0;
    //m_recChain = 0;
    //m_mcChain = 0;
    
    evt = 0;
    //rec = 0;
    //mc = 0;
    
    fileArr = 0;
    treeArr = 0;
    chainArr = 0;
}

#endif
