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
      sprintf(temp,"varCorr%1dp%1dv", iPlane, iView);
      m_chargeStrip[iPlane][iView] = new TGraphErrors(g_nDiv);
      m_chargeStrip[iPlane][iView]->SetName(temp);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	char name1[] = "tot0p0v0000fe";
	sprintf(name1,"tot%1dp%1dv%04dfe", iPlane, iView, iDiv);

	m_totHist[iPlane][iView][iDiv] = new TH1F(name1, name1, 200, 0, 200);

	char name2[] = "charge0p0v0000fe";
	sprintf(name2,"charge%1dp%1dv%04dfe", iPlane, iView, iDiv);

	m_chargeHist[iPlane][iView][iDiv] = new TH1F(name2, name2, 120, 0, 30);


      }

    }
  }
}

totCalib::~totCalib() 
{
  if(m_totFile == 0) return;

  m_totFile->cd();
  
  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {

    for(int iView = 0; iView != g_nView; ++iView) {

      m_totStrip[iPlane][iView]->Write(0, TObject::kOverwrite);
      m_chargeStrip[iPlane][iView]->Write(0, TObject::kOverwrite);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	m_totHist[iPlane][iView][iDiv]->Write(0, TObject::kOverwrite);
	m_chargeHist[iPlane][iView][iDiv]->Write(0, TObject::kOverwrite);

      }
    }
  }

  m_totFile->Close();
}


bool totCalib::setOutputFiles( const char* outputTxt, const char* outputXml, 
		     const char* outputRoot )
{
  m_txtOutput = outputTxt;
  m_xmlOutput = outputXml;//takuya

  m_totFile = new TFile(outputRoot, "RECREATE");
  if( m_totFile ) std::cout << "Open output root file: " << outputRoot << std::endl;
  else{
    std::cout << outputRoot << " can not be opened." << std::endl;
    return false;
  }
  return true;
}


int totCalib::setInputRootFiles(const char* digi, const char* recon)
{
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

  int nEvents, nRecon, nDigi;
  if(m_reconFile) {
    nRecon = (int) m_reconTree->GetEntries();
    cout << "No of events in " << recon << " : " << nRecon << endl;
    nEvents = nRecon;
  }
  if(m_digiFile) {
    nDigi = (int) m_digiTree->GetEntries();
    cout << "No of events in " << digi << " : " << nDigi << endl;
    nEvents = nDigi;
  }

  if(nDigi != nRecon) {
    std::cout << "No. of events in the digi file is not equal to no. of events in the recon file! abort!" << std::endl;
    exit(1);
  }

  return nEvents;
}

int totCalib::setInputRootFiles( TChain* digi, TChain* recon )
{

  m_reconTree = recon;
  m_reconEvent = 0;
  m_reconTree->SetBranchAddress("ReconEvent", &m_reconEvent);

  m_digiTree = digi;
  m_digiEvent = 0;
  m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvent);

  int nEvents, nRecon, nDigi;
  nRecon = (int) m_reconTree->GetEntries();
  cout << "No of events in recon tree: " << nRecon << endl;
  nEvents = nRecon;
  
  nDigi = (int) m_digiTree->GetEntries();
  cout << "No of events in digi tree: " << nDigi << endl;
  nEvents = nDigi;

  if(nDigi != nRecon) {
    std::cout << "No. of events in the digi file is not equal to no. of events in the recon file! abort!" << std::endl;
    exit(1);
  }

  return nEvents;
}


void totCalib::calibChargeScale( int nEvents )
{
  nEvents = 50000;

  analyzeEvent(nEvents);

  fitTot();
  fillXml();//takuya
}

void totCalib::analyzeEvent(int nEvents) 
{
  int mEvent = nEvents * 0.1;
  for(int iEvent = 0; iEvent != nEvents; ++iEvent) {
    
    if( iEvent >= mEvent ){
      std::cout << "# of events: " << iEvent << " " << nEvents << std::endl;
      mEvent += nEvents * 0.1;
    }

    m_reconTree->GetEntry(iEvent);
    m_digiTree->GetEntry(iEvent);

    assert(m_reconEvent != 0);
    assert(m_digiEvent != 0);

    if(! passCut()) continue;

    getTot();

    fillTot();
  }
  std::cout << "Data scan finished." << std::endl;

}

void totCalib::getTot()
{
  int noOfTkrDigis = m_digiEvent->getTkrDigiCol()->GetLast()+1;

  for(int i = 0; i != noOfTkrDigis; ++i) {
    const TkrDigi* tkrDigi = m_digiEvent->getTkrDigi(i);

    assert(tkrDigi != 0);

    int iLayer = tkrDigi->getBilayer();
    int iPlane = g_nLayer - 1 - iLayer;
    //std::cout << iLayer << ", " << iPlane <<  std::endl;

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

//  assert(nHitPlane == 6);

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

      int tot = findTot(planeId, viewId, iStrip);

      float charge = calcCharge(planeId, viewId, iStrip, tot);

      static int nStripPerGroup = g_nStrip / g_nDiv;

      m_totHist[planeId][viewId][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
      m_chargeHist[planeId][viewId][iStrip/nStripPerGroup]->Fill(charge*(-m_dir.z()));
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
  if( output )
    std::cout << "Open log file: " << m_txtOutput << std::endl;
  else
    std::cout << m_txtOutput << " cannot be opened." << std::endl;

  // define Gaussian convolved Laudau function.
  TF1 *ffit = new TF1( "langau", langaufun, 0, 30, 4 );
  std::cout << "Start fit." << std::endl;

  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {
    for(int iView = 0; iView != g_nView; ++iView) {
      cout << "Layer: " << iPlane << ", View: " << iView << ", FE:";
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv){
	cout << " " << iDiv;
	
	if(m_totHist[iPlane][iView][iDiv]->GetEntries() == 0) continue;

	// fit uncorrected tot for each strip
	float ave = m_totHist[iPlane][iView][iDiv]->GetMean();
	float rms = m_totHist[iPlane][iView][iDiv]->GetRMS();
	float area = m_totHist[iPlane][iView][iDiv]->GetEntries();
	ffit->SetRange( ave-2*rms, ave+3*rms );
	ffit->SetParameters( rms*0.2, ave*0.75, area, rms*0.4 );
	ffit->SetParNames( "Width", "MP", "Area", "GSigma" );
	m_totHist[iPlane][iView][iDiv]->Fit( "langau", "RQ" );

	//0:width(scale) 1:peak 2:total area 3:width(sigma)
	double *par = ffit->GetParameters();
	double *error = ffit->GetParErrors();

	float pos = float(iDiv);
	float errPos = 0.;

	float peak = float( *(par+1) );
	float errPeak = float( *(error+1) );

	float width = float( *(par+3) );
	float errWidth = float( *(error+3) );

	m_totStrip[iPlane][iView]->SetPoint(iDiv, pos, peak);
	m_totStrip[iPlane][iView]->SetPointError(iDiv, errPos, errPeak);

	output << "Uncorrected tot " << iPlane << ' ' << iView << ' ' << pos
	       << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	       << errWidth << endl;

	// fit charge for each strip
	ave = m_chargeHist[iPlane][iView][iDiv]->GetMean();
	rms = m_chargeHist[iPlane][iView][iDiv]->GetRMS();
	m_chargeHist[iPlane][iView][iDiv]->Fit("landau", "RQ", "", ave-2*rms, ave+3*rms);

	par = (m_chargeHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParameters();
	error = (m_chargeHist[iPlane][iView][iDiv]->GetFunction("landau"))->GetParErrors();

        peak = float( *(par+1) );
	errPeak = float( *(error+1) );

	width = float( *(par+3) );
	errWidth = float( *(error+3) );

	m_chargeStrip[iPlane][iView]->SetPoint(iDiv, pos, peak);
	m_chargeStrip[iPlane][iView]->SetPointError(iDiv, errPos, errPeak);

	if( peak != 0.0 )
	  m_chargeScale[iPlane][iView][iDiv] = 5.0 / peak;
	else
	  m_chargeScale[iPlane][iView][iDiv] = 1.0;
	
	output << "Charge " << iPlane << ' ' << iView << ' ' << pos
	       << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	       << errWidth << endl;
	
      }
      cout << endl;
    }
  }
}


bool totCalib::readTotConvFile(const char* dir, const char* runid)
{
  string filename;
  char fname[] = "/398000364/TkrTotGainNt_LayerY17_398000364.tnt";
  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {
    for(int iView = 0; iView != g_nView; ++iView) {
      filename = dir;
      if( iView == 0 )
	sprintf(fname,"/%s/TkrTotGainNt_LayerX%d_%s.tnt", runid, iPlane, runid);
      else
	sprintf(fname,"/%s/TkrTotGainNt_LayerY%d_%s.tnt", runid, iPlane, runid);
      filename += fname;
      if( !readTotConv( iPlane, iView, filename.c_str() ) )
	return false;
    }
  }
  return true;
}

bool totCalib::readTotConv(int layer, int view, const char* file)
{
  ifstream convFile(file);
  if(  !convFile ){
    std::cout << file << " cannot be opened." << std::endl;
    return false;
  }
  else std::cout << "Reading " << file << std::endl;
  for(int i = 0; i != 2; ++i) {
    string temp;
    getline(convFile, temp);
  }

  int stripId, feId;
  float gain, offset, quadra, chisq;
  bool display = false;

  while(convFile >> stripId >> feId >> offset >> gain >> quadra >> chisq) {
    if( display ){
      std::cout << stripId << " " << offset << " " << quadra << std::endl;
      display = false;
    }
    m_totOffset[layer][view][stripId] = offset;
    m_totGain[layer][view][stripId] = gain;
    m_totQuadra[layer][view][stripId] = quadra;
  }

  return true;
}

float totCalib::calcCharge(int planeId, TkrCluster::view viewId, int iStrip, 
			   int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  int layer = g_nLayer - planeId - 1;
  int view = (viewId == TkrCluster::X) ? 0 : 1;

  // TOT to charge conversion
  float p0 = m_totOffset[layer][view][iStrip];
  float p1 = m_totGain[layer][view][iStrip];
  float p2 = m_totQuadra[layer][view][iStrip];
  float charge = 0.0;
  if( p2 != 0.0){
    float quad = p1*p1 - 4*p2*(p0-time);
    if( quad > 0.0 )
      charge = ( -p1 + sqrt( quad ) ) / (2*p2);
  }
  else if( p1 != 0.0 )
    charge = ( time - p0 ) / p1;
  
  return charge;
}


void totCalib::fillXml()//takuya
{
  std::ofstream output(m_xmlOutput.c_str());
  output << "<!-- dtd for Tracker calibrations -->" << endl
	 << "<!ELEMENT tot (generic?, cuts?, tower*) >" << endl
	 <<"<!ELEMENT threshold (generic?, cuts?, tower*) >" << endl
	 <<"<!ELEMENT chargeScale (generic?, cuts?, tower*) >" << endl
	 <<"<!ELEMENT  generic  (inputSample) >"<<endl
	 <<"<!ATTLIST generic"<<endl
	 <<"   instrument  (LAT | BFEM | BTEM | EM | CU | TWR) #REQUIRED"<<endl
	 <<"   timestamp   NMTOKEN #IMPLIED"<<endl
	 <<"   runId       NMTOKEN #IMPLIED"<<endl
	 <<"   calType     NMTOKEN #IMPLIED"<<endl
	 <<"   DTDVersion  NMTOKEN 'v1r0'"<<endl
	 <<"   fmtVersion  NMTOKEN #IMPLIED"<<endl
	 <<"   creatorName NMTOKEN #IMPLIED"<<endl
	 <<"   creatorVersion NMTOKEN #IMPLIED "<<endl
	 <<" >"<<endl

	 <<"<!-- Description of events used as input to the procedure"<<endl
	 <<"     Start- and stop-times should be timestamps of earliest and"<<endl
	 <<"     latest events included in sample"<<endl
	 <<"-->"<<endl
	 <<"<!ELEMENT inputSample (#PCDATA) >"<<endl
	 <<"<!ATTLIST inputSample  startTime NMTOKEN #REQUIRED"<<endl
	 <<"                       stopTime  NMTOKEN #REQUIRED"<<endl
	 <<"                       triggers  NMTOKENS #REQUIRED"<<endl
	 <<"                       source    NMTOKENS #REQUIRED"<<endl
	 <<"                       mode      NMTOKEN  #REQUIRED>"<<endl
	 <<""<<endl
	 <<"<!ELEMENT cuts EMPTY>"<<endl
	 <<"<!ATTLIST cuts  tight       NMTOKEN #REQUIRED"<<endl
	 <<"                loose       NMTOKEN #REQUIRED"<<endl
	 <<"                expected    NMTOKEN #REQUIRED >"<<endl
	 <<""<<endl
	 <<""<<endl
	 <<"<!ELEMENT tower (uniplane+) >"<<endl
	 <<"<!ATTLIST tower"<<endl
	 <<"   row      NMTOKEN #REQUIRED"<<endl
	 <<"   col      NMTOKEN #REQUIRED"<<endl
	 <<"   hwserial NMTOKEN #IMPLIED >"<<endl
    
	 <<"<!ELEMENT uniplane ((stripTot+) | (stripScale+) | (stripThresh+) |"<<endl
	 <<"                    (gtfeScale+) | (gtfeThresh+) ) >"<<endl
	 <<"<!ATTLIST uniplane"<<endl
	 <<"  tray NMTOKEN #REQUIRED"<<endl
	 <<"  which (bot | top ) #REQUIRED >"<<endl<<endl

	 <<"<!ELEMENT stripTot EMPTY >"<<endl
	 <<"<!ATTLIST stripTot"<<endl
	 <<"   id        NMTOKEN #REQUIRED"<<endl
	 <<"   slope     NMTOKEN #REQUIRED"<<endl
	 <<"   intercept NMTOKEN #REQUIRED"<<endl
	 <<"   quad      NMTOKEN #REQUIRED"<<endl
	 <<"   chi2      NMTOKEN #REQUIRED"<<endl
	 <<"   df        NMTOKEN #REQUIRED >"<<endl


	 <<"<!ELEMENT stripScale EMPTY >"<<endl
	 <<"<!ATTLIST stripScale"<<endl
	 <<"   id          NMTOKEN #REQUIRED"<<endl
	 <<"   chargeScale NMTOKEN #IMPLIED >"<<endl
    

	 <<"<!ELEMENT stripThresh EMPTY >"<<endl
	 <<"<!ATTLIST stripThresh"<<endl
	 <<"   id        NMTOKEN #REQUIRED"<<endl
	 <<"   trg_thr      NMTOKEN #IMPLIED"<<endl
	 <<"   data_thr     NMTOKEN #IMPLIED>"<<endl


	 <<"<!ELEMENT gtfeScale EMPTY >"<<endl
	 <<"<!ATTLIST gtfeScale"<<endl
	 <<"  id           NMTOKEN #REQUIRED"<<endl
	 <<"  chargeScale  NMTOKEN #IMPLIED >"<<endl
    
	 <<"<!ELEMENT gtfeThresh EMPTY >"<<endl
	 <<"<!ATTLIST gtfeThresh"<<endl
	 <<"   id        NMTOKEN #REQUIRED"<<endl
	 <<"   trg_thr   NMTOKEN #IMPLIED"<<endl
	 <<"   data_thr  NMTOKEN #IMPLIED >"<<endl;


  output << "<TrackerCalib>" << endl
	 << " <tot>" << endl
	 << "  <tower row='!' col='!'>" << endl;
  for(int iPlane = 0; iPlane != g_nPlane; ++iPlane) {
    for(int iView = 0; iView != g_nView; ++iView) {
      string a;
      if(iView==0) a="top";
      else if(iView==1) a="bot";
      output   << "   <uniplane tray='" << iPlane << "' which='"
	       << a << "'>" << endl;
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	output << "    <gtfeScale id='" << iDiv << "' chargeScale='" 
               <<  m_chargeScale[iPlane][iView][iDiv] << "'/>" << endl;
      }
      output << "   </uniplane>" << endl; 
    }
  }
  output     << "  </tower>" << endl;  
  output     << " </tot>" << endl
	     << "</TrackerCalib>" << endl;
}


//-----------------------------------------------------------------------
//
//      Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//  Markus Friedl (Markus.Friedl@cern.ch)
//  Modified by Hiro Tajima (tajima@stanford.edu) for speed up.
//
//-----------------------------------------------------------------------

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 40.0;   // number of convolution steps
  Double_t sc =  2.0;   // convolution extends to +-sc Gaussian sigmas
  Double_t np2 = 20.0;  // number of convolution steps
  Double_t sc2 =  4.0;  // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland, fgauss;
  Double_t sum = 0.0;
  Double_t xm, dx;
  Double_t step;
  Double_t sigma = par[3];
  Double_t width = par[0];

  // MP shift correction
  mpc = par[1] - mpshift * width;
  
  // Convolution integral parameters
  xm = x[0];
  Double_t dxmax = sc * sigma;
  step = dxmax * 2 / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=0.5*step; dx<=dxmax; dx+=step) {
    fgauss = TMath::Gaus(0.0,dx,sigma);
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
  }
  sum *= step/width;
  
  Double_t dxmax2 = sc2 * sigma;
  step = (dxmax2-dxmax) * 2 / np2;
  Double_t sum2 = 0.0;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=dxmax+0.5*step; dx<=dxmax2; dx+=step) {
    fgauss = TMath::Gaus(0.0,dx,sigma);
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum2 += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum2 += fland * fgauss;
  }
  sum2 *= step/width;
  
  return (par[2] * (sum+sum2) * invsq2pi / sigma);
}
