#include <cmath>
#include <fstream>
#include <iostream>
#include "TF1.h"
#include "totCalib.h"

using std::string;
using std::cout;
using std::endl;

totCalib::totCalib(): m_reconFile(0), m_reconTree(0), 
    m_reconEvent(0), m_digiFile(0), m_digiTree(0),
    m_digiEvent(0), m_totFile(0)
{
  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {

    for(int iView = 0; iView != g_nView; ++iView) {
      char name[] = "var00";
      sprintf(name,"var%1d%1d", iPlane, iView);
      m_totStrip[iPlane][iView] = new TGraphErrors(g_nDiv);
      m_totStrip[iPlane][iView]->SetName(name);

      char temp[] = "varCorr00";
      sprintf(temp,"varCorr%1d%1d", iPlane, iView);
      m_totCorrStrip[iPlane][iView] = new TGraphErrors(g_nDiv);
      m_totCorrStrip[iPlane][iView]->SetName(temp);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	char name1[] = "tot000000";
	sprintf(name1,"tot%1d%1d%04d", iPlane, iView, iDiv);

	m_totHist[iPlane][iView][iDiv] = new TH1F(name1, name1, 250, 0, 500);

	char name2[] = "totCorr000000";
	sprintf(name2,"totCorr%1d%1d%04d", iPlane, iView, iDiv);

	m_totCorrHist[iPlane][iView][iDiv] = new TH1F(name2, name2, 120, 0, 30);


      }

    }
  }

  readTotCorr(1, 0, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_x1.txt");
  readTotCorr(1, 1, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_y1.txt");
  readTotCorr(2, 0, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_x2.txt");
  readTotCorr(2, 1, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_y2.txt");
  readTotCorr(3, 0, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_x3.txt");
  readTotCorr(3, 1, "/nfs/farm/g/glast/u03/EM2003/rootFiles/em_v1r030302p5/tot//chargeInjection_y3.txt");

}

totCalib::~totCalib() 
{
  if(m_totFile == 0) return;

  m_totFile->cd();
  
  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {

    for(int iView = 0; iView != g_nView; ++iView) {

      m_totStrip[iPlane][iView]->Write(0, TObject::kOverwrite);
      m_totCorrStrip[iPlane][iView]->Write(0, TObject::kOverwrite);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	m_totHist[iPlane][iView][iDiv]->Write(0, TObject::kOverwrite);
	m_totCorrHist[iPlane][iView][iDiv]->Write(0, TObject::kOverwrite);

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
    std::cout << "No. of events in the digi file is not equal to no. of events in the recon file! abort!" << std::endl;
    exit(1);
  }

  //  nEvent = 10000;

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
    std::cout << "No. of events in the digi file is not equal to no. of events in the recon file! abort!" << std::endl;
    exit(1);
  }

  //  nEvent = 10000;

  analyzeEvent(nEvent);

  fitTot();
}

void totCalib::analyzeEvent(int nEvent) 
{
  for(int iEvent = 0; iEvent != nEvent; ++iEvent) {

    m_reconTree->GetEntry(iEvent);
    m_digiTree->GetEntry(iEvent);

    assert(m_reconEvent != 0);
    assert(m_digiEvent != 0);

    if(! passCut()) continue;

    getTot();

    fillTot();
  }
}

void totCalib::getTot()
{
  int noOfTkrDigis = m_digiEvent->getTkrDigiCol()->GetLast()+1;

  for(int i = 0; i != noOfTkrDigis; ++i) {
    const TkrDigi* tkrDigi = m_digiEvent->getTkrDigi(i);

    assert(tkrDigi != 0);

    int iLayer = tkrDigi->getBilayer();
    int iPlane = g_nLayer - 1 - iLayer;

    GlastAxis::axis view = tkrDigi->getView();
    if(view == GlastAxis::X) {
      m_totX[iPlane][0] = tkrDigi->getToT(0);
      m_totX[iPlane][1] = tkrDigi->getToT(1);
    }
    else if(view == GlastAxis::Y) {
      m_totY[iPlane][0] = tkrDigi->getToT(0);
      m_totY[iPlane][1] = tkrDigi->getToT(1);
    }
    else {
      cout << "Unknown view: " << int(view) << endl;
    }
  } 
}
 
void totCalib::retrieveCluster()
{
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
 
  assert(tkrRecon != 0);

  TObjArray* siClusterCol = tkrRecon->getClusterCol();

  int noOfTkrClusters = siClusterCol->GetLast()+1;
  for(int i = 0; i != noOfTkrClusters; ++i) {
    TkrCluster* cluster = dynamic_cast<TkrCluster*>(siClusterCol->At(i));
    m_cluster[cluster->getId()] = cluster;
  }
}

int totCalib::findTot(int planeId, TkrCluster::view viewId, int stripId)
{
  if(planeId != 2) {
    if(stripId < g_nStrip/2) {
      if(viewId == TkrCluster::X) {
	return m_totX[planeId][0];
      }
      else {
	return m_totY[planeId][0];
      }
    }
    else {
      if(viewId == TkrCluster::X) {
	return m_totX[planeId][1];
      }
      else {
	return m_totY[planeId][1];
      }
    }
  }
  else {
    if(stripId < g_nStrip/g_nFecd * 4) {
      if(viewId == TkrCluster::X) {
	return m_totX[planeId][0];
      }
      else {
	return m_totY[planeId][0];
      }
    }
    else {
      if(viewId == TkrCluster::X) {
	return m_totX[planeId][1];
      }
      else {
	return m_totY[planeId][1];
      }
    }
  }
}

void totCalib::fillTot() 
{
  retrieveCluster();

  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();

  TObjArray* tracks = tkrRecon->getTrackCol();
  TkrKalFitTrack* tkrTrack = dynamic_cast<TkrKalFitTrack*>(tracks->At(0));

  int nHitPlane = tkrTrack->getNumHits();

  assert(nHitPlane == 6);

  for(int iPlane = 0; iPlane != nHitPlane; ++iPlane) {

    const TkrHitPlane* plane = tkrTrack->getHitPlane(iPlane);

    std::map<int, TkrCluster*>::const_iterator itr = m_cluster.find(plane->getIdHit());

    assert(itr != m_cluster.end());

    TkrCluster* cluster = itr->second;

    int planeId = cluster->getPlane();
    TkrCluster::view viewId = cluster->getView();

    // require only a single strip
    if(cluster->getSize() != 1) continue;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip) {

      float tot = findTot(planeId, viewId, iStrip);

      float totCorr = correctTot(planeId, viewId, iStrip, tot);

      static int nStripPerGroup = g_nStrip / g_nDiv;

      m_totHist[planeId][viewId][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
      m_totCorrHist[planeId][viewId][iStrip/nStripPerGroup]->Fill(totCorr*(-m_dir.z()));
    }
  }
}

bool totCalib::passCut() 
{
    TkrRecon* tkrRecon = m_reconEvent->getTkrRecon(); 
    assert(tkrRecon != 0);

    TObjArray* vertices = tkrRecon->getVertexCol();

    // select only 1 track event
    if(vertices->GetLast()+1 != 1) return false;

    TkrVertex* tkrVertex = dynamic_cast<TkrVertex*>(vertices->At(0));
    if(tkrVertex) {
      m_pos = tkrVertex->getPosition();
      m_dir = tkrVertex->getDirection();

      if(m_dir.Z() > -0.9) return false;
    }

    return true;
}

void totCalib::fitTot()
{  
  std::ofstream output(m_txtOutput.c_str());

  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {
    for(int iView = 0; iView != g_nView; ++iView) {
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	
	if(m_totHist[iPlane][iView][iDiv]->GetEntries() == 0) continue;

	// fit uncorrected tot for each strip
	float ave = m_totHist[iPlane][iView][iDiv]->GetMean();
	float rms = m_totHist[iPlane][iView][iDiv]->GetRMS();
	m_totHist[iPlane][iView][iDiv]->Fit("landau", "", "", ave-2*rms, ave+3*rms);

	double* par = (m_totHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParameters();

	double* error = (m_totHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParErrors();

	float pos = (iDiv);
	float errPos = 0.;

	float peak = float( *(par+1) );
	float errPeak = float( *(error+1) );

	float width = float( *(par+2) );
	float errWidth = float( *(error+2) );

	m_totStrip[iPlane][iView]->SetPoint(iDiv, pos, peak);
	m_totStrip[iPlane][iView]->SetPointError(iDiv, errPos, errPeak);

	output << "Uncorrected tot " << iPlane << ' ' << iView << ' ' << pos
	       << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	       << errWidth << endl;

	// fit corrected tot for each strip
	ave = m_totCorrHist[iPlane][iView][iDiv]->GetMean();
	rms = m_totCorrHist[iPlane][iView][iDiv]->GetRMS();
	m_totCorrHist[iPlane][iView][iDiv]->Fit("landau", "", "", ave-2*rms, ave+3*rms);

	par = (m_totCorrHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParameters();

	error = (m_totCorrHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParErrors();

	pos = (iDiv);
	errPos = 0.;

        peak = float( *(par+1) );
	errPeak = float( *(error+1) );

	width = float( *(par+2) );
	errWidth = float( *(error+2) );

	m_totCorrStrip[iPlane][iView]->SetPoint(iDiv, pos, peak);
	m_totCorrStrip[iPlane][iView]->SetPointError(iDiv, errPos, errPeak);

	output << "Corrected tot " << iPlane << ' ' << iView << ' ' << pos
	       << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	       << errWidth << endl;
	
      }

    }
  }
}

void totCalib::readTotCorr(int layer, int view, const char* file)
{
  ifstream corrFile(file);
  for(int i = 0; i != 14; ++i) {
    string temp;
    getline(corrFile, temp);
  }

  int stripId;
  float gain, offset;

  while(corrFile >> stripId >> gain >> offset) {
    m_totGain[layer][view][stripId] = gain;
    m_totOffset[layer][view][stripId] = offset;
  }
}

float totCalib::correctTot(int planeId, TkrCluster::view viewId, int iStrip, 
			   int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  int layer = g_nLayer - planeId - 1;
  int view = (viewId == TkrCluster::X) ? 0 : 1;

  // TOT correction
  return (time  - m_totOffset[layer][view][iStrip]) / m_totGain[layer][view][iStrip];
}
