#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <iomanip>
#include "TF1.h"
#include "totCalib.h"

#include "commonRootData/idents/TowerId.h"
#include "digiRootData/TkrDigi.h"

using std::string;
using std::cout;
using std::endl;

// to put all the data into one tower, for testing
const bool overlayTowers = false; 

// Some special histograms
TH1F *outsideTowers, *outsidePlanes;

totCalib::totCalib(): m_reconFile(0), m_reconTree(0), 
m_reconEvent(0), m_digiFile(0), m_digiTree(0),
m_digiEvent(0), m_totFile(0)
{
    // set up the vectors
    m_totStrip.setDims(g_nTower, g_nLayer, g_nView);
    m_totHist.setDims(g_nTower, g_nLayer, g_nView, g_nDiv);

    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111);

    m_peaks = new TH1F("FittedPeaks", "Fitted Peaks", 90, 0, 3.);
    m_widths = new TH1F("FittedWidths", "Fitted Widths", 40, 0, 1.);
    m_chisqs = new TH1F("ChiSquares", "ChiSquares of Fits", 100, 0., 5);

    outsideTowers = new TH1F("outsideTowers", "hits in Towers outside current config", 18, -1.5, 16.5);
    outsidePlanes = new TH1F("outsidePlanes", "hits in planes outside current config", 38, -1.5, 36.5);


    for(int iTower = 0; iTower<g_nTower; ++iTower) {
        for(int iLayer = 0; iLayer != g_nLayer; ++iLayer) {
            for(int iView = 0; iView != g_nView; ++iView) {
                /*
                char name[] = "var00:00-0";
                sprintf(name,"var%2d:%2d-%1d", iTower, iLayer, iView);
                int index = m_totStrip.getIndex(iTower, iLayer,iView);
                m_totStrip[index] = new TGraphErrors(g_nDiv);
                m_totStrip[index]->SetName(name);
                */

                for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
                    char name1[] = "tot00:00-0-0000";
                    int index1 = m_totHist.getIndex(iTower, iLayer, iView, iDiv);
                    sprintf(name1,"tot%2d-%2d-%1d-%04d", iTower, iLayer, iView, iDiv);

                    m_totHist[index1] = new TH1F(name1, name1, 50, 0, 5.);
                }
            }
        }
    }
    std::cout << "Histograms initialized" << std::endl;
}

totCalib::~totCalib() 
{
    if(m_totFile == 0) return;
    m_totFile->cd();
    outsideTowers->Write(0, TObject::kOverwrite);
    outsidePlanes->Write(0, TObject::kOverwrite);

    m_peaks->Write(0, TObject::kOverwrite);
    m_widths->Write(0, TObject::kOverwrite);
    m_chisqs->Write(0, TObject::kOverwrite);

    for(int iTower = 0; iTower<g_nTower; ++iTower) {
        for(int iLayer = 0; iLayer != g_nLayer; ++iLayer) {
            for(int iView = 0; iView != g_nView; ++iView) {
                int index = m_totStrip.getIndex(iTower, iLayer, iView);
                if (m_totStrip[index]) m_totStrip[index]->Write(0, TObject::kOverwrite);

                for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
                    int index1 = m_totHist.getIndex(iTower, iLayer, iView, iDiv);
                    m_totHist[index1]->Write(0, TObject::kOverwrite);
                }
            }
        }
    }
    m_totFile->Close();
}

void totCalib::genTot(const char* digi, const char* recon, 
                      const char* outputTxt, const char* outputRoot)
{
    m_txtOutput = outputTxt;
    m_totFile = new TFile(outputRoot, "RECREATE");

    m_reconFile = new TFile(recon, "READ");
    if(m_reconFile->IsZombie()) {
        m_reconFile = 0;
        std::cout << "recon file " << recon << " does not exist! abort!" <<
            std::endl;
        exit(1);
    }

    if(m_reconFile) {
        m_reconTree = (TTree*) m_reconFile->Get("Recon");
        m_reconEvent = 0;
        m_reconTree->SetBranchAddress("ReconEvent", &m_reconEvent);
    }

    m_digiFile = new TFile(digi, "READ");
    if(m_digiFile->IsZombie()) {
        m_digiFile = 0;
        std::cout << "digi file " << digi << " does not exist! abort!" <<
            std::endl;
        exit(1);
    }

    if(m_digiFile) {
        m_digiTree = (TTree*) m_digiFile->Get("Digi");
        m_digiEvent = 0;
        m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvent);
    }

    int nEvent, nRecon, nDigi;
    if(m_reconFile) {
        nRecon = (int) m_reconTree->GetEntries();
        cout << "No of events in " << recon << " : " << nRecon << endl;
        nEvent = nRecon;
    }
    if(m_digiFile) {
        nDigi = (int) m_digiTree->GetEntries();
        cout << "No of events in " << digi << " : " << nDigi << endl;
        nEvent = nDigi;
    }

    if(nDigi != nRecon) {
        std::cout << "No. of events in the digi file is not equal to "
            << "no. of events in the recon file! abort!" << std::endl;
        exit(1);
    }

    analyzeEvent(nEvent);

    fitTot();
}

void totCalib::genTot(TChain* digi, TChain* recon, 
                      const char* outputTxt, const char* outputRoot)
{
    m_txtOutput = outputTxt;
    m_totFile = new TFile(outputRoot, "RECREATE");

    m_reconTree = recon;
    m_reconEvent = 0;
    m_reconTree->SetBranchAddress("ReconEvent", &m_reconEvent);

    m_digiTree = digi;
    m_digiEvent = 0;
    m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvent);

    int nEvent, nRecon, nDigi;
    nRecon = (int) m_reconTree->GetEntries();
    cout << "No of events in recon tree: " << nRecon << endl;
    nEvent = nRecon;

    nDigi = (int) m_digiTree->GetEntries();
    cout << "No of events in digi tree: " << nDigi << endl;
    nEvent = nDigi;

    if(nDigi != nRecon) {
        std::cout << "No. of events in the digi file is not equal to"
            << " no. of events in the recon file! abort!" << std::endl;
        exit(1);
    }

    //  nEvent = 10000;

    analyzeEvent(nEvent);

    fitTot();
}

void totCalib::analyzeEvent(int nEvent) 
{
    for(int iEvent = 0; iEvent != nEvent; ++iEvent) {
        if(iEvent%1000==0) std::cout << "Processing event " << iEvent << std::endl;

        m_reconTree->GetEntry(iEvent);
        m_digiTree->GetEntry(iEvent);

        assert(m_reconEvent != 0);
        assert(m_digiEvent != 0);

        TkrRecon* tkrRecon = m_reconEvent->getTkrRecon(); 
        assert(tkrRecon != 0);

        if(! passCut(tkrRecon)) continue;

        fillTot(tkrRecon);
    }
}


void totCalib::fillTot(TkrRecon* tkrRecon) 
{
    TObjArray* clusterCol = tkrRecon->getClusterCol();

    TObjArray* tracks = tkrRecon->getTrackCol();
    TkrTrack* track = dynamic_cast<TkrTrack*>(tracks->At(0));

    int size = track->GetEntriesFast();

    int nStripPerGroup = g_nStrip / g_nDiv;

    double costh = -m_dir.Z();

    int i;
    for (i=0; i<size; ++i) {

        TkrTrackHit* iter = dynamic_cast<TkrTrackHit*>(track->At(i));
        const TkrCluster* cluster = ((*iter).getClusterPtr());
        if (!cluster) continue;

        int layer = cluster->getLayer();
        //TkrCluster::view viewId = cluster->getView();

        // require only a single strip
        if(cluster->getSize() != 1) continue;
        commonRootData::TkrId id = cluster->getTkrId();
        int plane  = cluster->getPlane();
        int view   = id.getView();
        int towerX = id.getTowerX();
        int towerY = id.getTowerY();
        TowerId towerId;
        int tower = TowerId(towerX, towerY).id();
        //int tower = towerX + 4*towerY;
        if (overlayTowers) tower = 0;
        if( tower>=16 || tower < 0) {
            int truncTower = std::max(-1, std::min(16, tower));
            outsideTowers->Fill(truncTower);
            continue;
        }

        if (plane>=g_nPlane || plane < 0) {
            int truncPlane = std::max(-1, std::min(36, plane));
            outsidePlanes->Fill(truncPlane);
            continue;
        } 
        float tot = cluster->getMips();
        int   end = cluster->getEnd();
        if (end>1) continue;

        int iStrip = cluster->getFirstStrip();
        int fStrip = cluster->getLastStrip();
        for( ; iStrip != fStrip+1; ++iStrip) {
            m_totHist(tower, layer, view, iStrip/nStripPerGroup)->Fill(tot*costh); 
        }
    }
}


bool totCalib::passCut(TkrRecon* tkrRecon) 
{

    // should be tracks, not vertices
    TObjArray* tracks = tkrRecon->getTrackCol();

    // select only 1 track event
    if(tracks->GetEntriesFast() != 1) return false;

    // reject really short tracks?
    //if (tkrTrack->getNumHits()<8) return false;
    
    TkrTrack* track = dynamic_cast<TkrTrack*>(tracks->At(0));
    if(track) {
        m_pos = track->getInitialPosition();
        m_dir = track->getInitialDirection();

        if(m_dir.Z() > -0.9) return false;
    }
  
    return true;
}

void totCalib::fitTot()
{  
    std::ofstream output(m_txtOutput.c_str());
    int nFit = 0;
    for (int i=0; i<(int)m_totStrip.size(); ++i) {
        if(m_totHist[i]->GetEntries()>0) ++nFit;
    }

    std::cout << "Fitting " << nFit << " histograms " << std::endl;
    for(int iTower = 0; iTower < g_nTower; ++iTower) {
        for(int iLayer = 0; iLayer != g_nLayer; ++iLayer) {
            for(int iView = 0; iView != g_nView; ++iView) {
                int index = m_totStrip.getIndex(iTower, iLayer, iView);
                std::vector<float> x, y, xerr, yerr;
                int nPoints = 0;
                for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
                    int index1 = m_totHist.getIndex(iTower, iLayer, iView, iDiv);
                    if(m_totHist[index1]->GetEntries() == 0) continue;

                    // fit corrected tot for each strip
                    float ave = m_totHist[index1]->GetMean();
                    float rms = m_totHist[index1]->GetRMS();
                    m_totHist[index1]->Fit("landau", "Q", "", ave-2*rms, ave+3*rms);

                    double* par = (m_totHist[index1]->GetFunction("landau"))->GetParameters();
                    double* error = (m_totHist[index1]->GetFunction("landau"))->GetParErrors();
                    double  chisq = (m_totHist[index1]->GetFunction("landau"))->GetChisquare();
                    int     NDoF  = (m_totHist[index1]->GetFunction("landau"))->GetNDF();


                    double pos      = iDiv;
                    double errPos   = 0.;
                    double peak     = par[1];
                    double errPeak  = error[1];
                    double width    = par[2];
                    double errWidth = error[2];

                    x.push_back(iDiv);
                    y.push_back(peak);
                    xerr.push_back(0);
                    yerr.push_back(errPeak);
                    nPoints++;

                    /*
                    m_totStrip[index]->SetPoint(iDiv, pos, peak);
                    m_totStrip[index]->SetPointError(iDiv, errPos, errPeak);
                    */
                    m_peaks->Fill(peak);
                    m_widths->Fill(width);
                    m_chisqs->Fill(chisq/std::max(1, NDoF));

                    output << "tot " << iTower << ": " << iLayer << ' ' << iView << ' ' << pos << '\t' 
                        << std::setw(10) << std::setprecision(4) << peak 
                        << std::setw(9)  << std::setprecision(3) << errPeak 
                        << std::setw(12) << std::setprecision(4) << width 
                        << std::setw(9)  << std::setprecision(3) << errWidth 
                        << std::setw(11) << std::setprecision(4) << chisq/std::max(1, NDoF) 
                        << endl;
                }
                if (nPoints==0) continue;
                char name[] = "var00:00-0";
                sprintf(name,"var%2d:%2d-%1d", iTower, iLayer, iView);
                //int index = m_totStrip.getIndex(iTower, iLayer,iView);
                m_totStrip[index] = new TGraphErrors(nPoints, &x[0], &y[0], &xerr[0], &yerr[0]);
                m_totStrip[index]->SetName(name);

            }
        }
    }
}
