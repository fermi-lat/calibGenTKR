
#include <cmath>
#include <ctime>

#include "totCalib.h"
 
using std::string;
using std::cout;
using std::endl;

XERCES_CPP_NAMESPACE_USE


totCalib::totCalib(): m_reconFile(0), m_reconTree(0), 
    m_reconEvent(0), m_digiFile(0), m_digiTree(0),
    m_digiEvent(0), m_totFile(0)
{

  m_tower_row=0;
  m_tower_col=0;
  m_tower_serial="TKrFMX";
  m_version="1.9";
  m_first_run = 999999999;
  m_last_run = 0;
  m_tot_runid = "0";
  m_startTime="01/20 2005, 22:52 GMT";
  m_stopTime="01/20 2005, 22:52 GMT";

  // get time 
  time_t rawtime;
  time( &rawtime );
  // format time into a null-terminated string
  char nts[] = "050124-000000";
  size_t ntsmax;
  strftime( nts, ntsmax, "%y%m%d", gmtime( &rawtime ) );
  m_dateStamp = nts;
  strftime( nts, ntsmax, "%H%M%S", gmtime( &rawtime ) );
  m_timeStamp = nts;

  for(int layer = 0; layer != g_nPlane; ++layer) {

    for(int iView = 0; iView != g_nView; ++iView) {
      char name[] = "var00";
      sprintf(name,"var%1d%1d", layer, iView);
      m_totStrip[layer][iView] = new TGraphErrors(g_nDiv);
      m_totStrip[layer][iView]->SetName(name);

      char temp[] = "varCorr00";
      sprintf(temp,"varCorr%1dp%1dv", layer, iView);
      m_chargeStrip[layer][iView] = new TGraphErrors(g_nDiv);
      m_chargeStrip[layer][iView]->SetName(temp);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	char name1[] = "tot0p0v0000fe";
	sprintf(name1,"tot%1dp%1dv%04dfe", layer, iView, iDiv);

	m_totHist[layer][iView][iDiv] = new TH1F(name1, name1, 100, 0, 200);

	char name2[] = "charge0p0v0000fe";
	sprintf(name2,"charge%1dp%1dv%04dfe", layer, iView, iDiv);

	m_chargeHist[layer][iView][iDiv] = new TH1F(name2, name2, 200, 0, 20);


      }

    }
  }
}

totCalib::~totCalib() 
{
  if(m_totFile == 0) return;

  m_totFile->cd();
  
  for(int layer = 0; layer != g_nPlane; ++layer) {

    for(int iView = 0; iView != g_nView; ++iView) {

      m_totStrip[layer][iView]->Write(0, TObject::kOverwrite);
      m_chargeStrip[layer][iView]->Write(0, TObject::kOverwrite);

      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {

	m_totHist[layer][iView][iDiv]->Write(0, TObject::kOverwrite);
	m_chargeHist[layer][iView][iDiv]->Write(0, TObject::kOverwrite);

      }
    }
  }

  m_totFile->Close();
}


bool totCalib::setOutputFiles( const char* outputDir )
{
  m_outputDir = outputDir;

  std::string filename;
  char fname[] = "/TE603_050121-0000.root";

  filename = m_outputDir;
  sprintf( fname, "/TE603_%s-%s.log", m_dateStamp.c_str(), m_timeStamp.c_str() );
  filename += fname;
  m_log.open( filename.c_str() );
  if( m_log )
    std::cout << "Open log file: " << filename << std::endl;
  else{
    std::cout << filename << " cannot be opened." << std::endl;
    return false;
  }

  filename = m_outputDir;
  sprintf( fname, "/TE603_%s-%s.root", m_dateStamp.c_str(), m_timeStamp.c_str() );
  filename += fname;
  m_totFile = new TFile( filename.c_str(), "RECREATE" );
  if( m_totFile ){
    std::cout << "Open output root file: " << filename << std::endl;
    m_log << "Output root file: " << filename << std::endl;
  }
  else{
    std::cout << filename << " can not be opened." << std::endl;
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


int totCalib::setInputRootFiles( const char* rootDir, const char* reconDir, 
				 const std::vector<std::string>& runIds )
{

  std::string digiFile, reconFile;
  char fname[] = "/398000362/calib-v1r0/grRoot/recon-EM2-v1r0_398000362_recon_RECON.root";
  TChain* digiChain = new TChain("Digi");
  TChain* reconChain = new TChain("Recon");

  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {

    int runid = atoi( (*run).c_str() );
    if( runid < m_first_run ) m_first_run = runid;
    if( runid > m_last_run ) m_last_run = runid;

    digiFile = rootDir;
    sprintf(fname,"/%d/grRoot/digitization-EM2-v1r0_%d_digi_DIGI.root",
	    runid, runid);
    digiFile += fname;
    std::cout << "open digi file: " << digiFile << endl;
    m_log << "digi file: " << digiFile << endl;
    digiChain->Add( digiFile.c_str() );

    reconFile = rootDir;
    sprintf(fname,"/%d/%s/grRoot/recon-EM2-v1r0_%d_recon_RECON.root",
	    runid, reconDir, runid);
    reconFile += fname;
    std::cout << "open recon file: " << reconFile << endl;
    m_log << "recon file: " << reconFile << endl;
    reconChain->Add( reconFile.c_str() );

  }
  return setInputRootFiles( digiChain, reconChain );

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


bool totCalib::readRcReports( const char* reportDir, 
				 const std::vector<std::string>& runIds )
{
  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {
    std::string reportFile = reportDir;
    reportFile += "/" + *run;
    reportFile += "/rcReport.out";
    std::cout << "open rcReport file: " << reportFile << endl;
    m_log << "rcReport file: " << reportFile << endl;
    if( !parseRcReport( reportFile.c_str() ) ) return false;
  }
  return true;
}


bool totCalib::parseRcReport( const char* reportFile )
{
  xml::XmlParser* parsercReport = new xml::XmlParser(true);
  DOMDocument* docrcReport = 0;
  try{
    docrcReport = parsercReport -> parse(reportFile);
  }
  catch (xml::ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    delete parsercReport;
    return false;
  }

  if (docrcReport != 0){//successful
    std::cout <<  reportFile << " is successfully parsed" << std::endl;
    
    //look up attributes
    string timeStamp;
    DOMElement* rcElt = docrcReport -> getDocumentElement();
    try {
      timeStamp = xml::Dom::getAttribute(rcElt, "timestamp");
    }
    catch (xml::DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    std::vector<std::string> keywords, values;
    keywords.push_back("RunId");
    keywords.push_back("StartTime");
    keywords.push_back("EndTime");

    for( int i=0; i<keywords.size(); i++){
      DOMElement* childElt 
	= xml::Dom::findFirstChildByName( rcElt, keywords[i].c_str() );
      try {
	values.push_back( xml::Dom::getTextContent(childElt) );
      }
      catch (xml::DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
    }
    int runid = atoi( values[0].c_str() );
    if( runid != m_first_run && runid != m_last_run ) return true;
    if( runid == m_first_run ) getDate( values[1].c_str(), m_startTime );
    if( runid == m_last_run ) getDate( values[2].c_str(), m_stopTime );
    
    return true;
  }
  return false;

}


void totCalib::getDate( const char* str, std::string& sdate )
{
  std::vector<std::string> strings;  
  std::string strs, schar;

  for(int i=0; str[i]; i++){
    if( str[i] == '(' || str[i] == ' ' ) continue;
    if( str[i] == ')' || str[i] == ',' ){
      strings.push_back( strs );
      strs.erase();
    }
    else{
      schar = str[i];
      strs.insert( strs.size(), schar );
    }
  }

  sdate = strings[1]; // month
  sdate += "/" + strings[2]; // day
  sdate += " " + strings[0]; // year
  sdate += ", " + strings[3]; // hour
  if( strings[4].size() == 1 )
    sdate += ":0" + strings[4]; // minutes
  else
    sdate += ":" + strings[4]; // minutes
  sdate += " GMT";
  std::cout << "date: " << sdate << std::endl;

}

void totCalib::calibChargeScale( int nEvents )
{
  //nEvents = 50000;

  analyzeEvent(nEvents);

  fitTot();
  fillXml();//takuya
}

void totCalib::analyzeEvent(int nEvents) 
{
  int mEvent = nEvents * 0.01;
  for(int iEvent = 0; iEvent != nEvents; ++iEvent) {
    
    if( iEvent >= mEvent ){
      std::cout << "# of events: " << iEvent << " " << nEvents 
		<< ",  " << iEvent*101/nEvents << "%" << std::endl;
      if( mEvent > nEvents*0.095 ) mEvent += nEvents * 0.1;
      else mEvent += nEvents * 0.01;
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

    int layer = g_nLayer - planeId - 1;
    int view = (viewId == TkrCluster::X) ? 0 : 1;

    // require only a single strip
    if(cluster->getSize() != 1) continue;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip) {

      int tot = findTot(planeId, viewId, iStrip);
      if( tot == 0 ) continue;

      float charge = calcCharge(layer, view, iStrip, tot);

      static int nStripPerGroup = g_nStrip / g_nDiv;

      m_totHist[layer][view][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
      m_chargeHist[layer][view][iStrip/nStripPerGroup]->Fill(charge*(-m_dir.z()));
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
  // define Gaussian convolved Laudau function.
  TF1 *ffit = new TF1( "langau", langaufun, 0, 30, 4 );
  ffit->SetParNames( "Width", "MP", "Area", "GSigma" );
  std::cout << "Start fit." << std::endl;

  for(int layer = 0; layer != g_nPlane; ++layer) {
    for(int iView = 0; iView != g_nView; ++iView) {
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv){
	std::cout << "Layer: " << layer << ", View: " << iView 
		  << ", FE: " << iDiv << std::endl;
	m_chargeScale[layer][iView][iDiv] = 1.0;

	// fit uncorrected tot for each strip
	float area = m_totHist[layer][iView][iDiv]->GetEntries();
	float ave = m_totHist[layer][iView][iDiv]->GetMean();
	float rms = m_totHist[layer][iView][iDiv]->GetRMS();
	//std::cout << area << " " << ave << " " << rms << std::endl;
	if( area<100 || ave==0.0 || rms==0.0 ){ 
	  std::cout << "Layer: " << layer << ", View: " << iView 
		    << ", FE: " << iDiv << ", Entries: " << area
		    << ", Mean: " << ave << ", RMS: " << rms 
		    << " skipped." << std::endl;
	  m_log << "Layer: " << layer << ", View: " << iView 
		    << ", FE: " << iDiv << ", Entries: " << area
		    << ", Mean: " << ave << ", RMS: " << rms 
		    << " skipped." << std::endl;
	  continue;
	}

	ffit->SetParLimits( 0, 0.0, rms );
	ffit->SetParLimits( 1, 0.0, ave*2 );
	ffit->SetParLimits( 2, 0.0, area*0.4 );
	ffit->SetParLimits( 3, 0.0, rms );
	ffit->SetRange( ave-2*rms, ave+3*rms );
	ffit->SetParameters( rms*0.2, ave*0.75, area*0.1, rms*0.4 );
	//m_totHist[layer][iView][iDiv]->Fit( "langau", "RBQ" );

	//0:width(scale) 1:peak 2:total area 3:width(sigma)
	Double_t *par = ffit->GetParameters();
	Double_t *error = ffit->GetParErrors();

	float pos = float(iDiv);
	float errPos = 0.;

	float peak = float( *(par+1) );
	float errPeak = float( *(error+1) );

	float width = float( *(par+3) );
	float errWidth = float( *(error+3) );

	m_totStrip[layer][iView]->SetPoint(iDiv, pos, peak);
	m_totStrip[layer][iView]->SetPointError(iDiv, errPos, errPeak);

	/*m_log << "Uncorrected tot " << layer << ' ' << iView << ' ' << pos
	       << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	       << errWidth << endl; */

	// fit charge for each strip
	area = m_chargeHist[layer][iView][iDiv]->GetEntries();
	ave = m_chargeHist[layer][iView][iDiv]->GetMean();
	rms = m_chargeHist[layer][iView][iDiv]->GetRMS();
	//std::cout << area << " " << ave << " " << rms << std::endl;
	if( area<100 || ave==0.0 || rms==0.0 ){ 
	  m_log << "Layer: " << layer << ", View: " << iView 
		    << ", FE: " << iDiv << ", Entries: " << area
		    << ", Mean: " << ave << ", RMS: " << rms 
		    << " skipped." << std::endl;
	  continue;
	}

	ffit->SetParLimits( 0, 0.0, rms );
	ffit->SetParLimits( 1, 0.0, ave*2 );
	ffit->SetParLimits( 2, 0.0, area*0.4 );
	ffit->SetParLimits( 3, 0.0, rms );
	ffit->SetRange( ave-2*rms, ave+3*rms );
	ffit->SetParameters( rms*0.2, ave*0.75, area*0.1, rms*0.4 );
	m_chargeHist[layer][iView][iDiv]->Fit( "langau", "RBQ" );

	//0:width(scale) 1:peak 2:total area 3:width(sigma)
	par = ffit->GetParameters();
	error = ffit->GetParErrors();
	//par = (m_chargeHist[layer][iView][iDiv]->GetFunction("landau"))->GetParameters();
	//error = (m_chargeHist[layer][iView][iDiv]->GetFunction("landau"))->GetParErrors();

        peak = float( *(par+1) );
	errPeak = float( *(error+1) );

	width = float( *(par+3) );
	errWidth = float( *(error+3) );

	m_chargeStrip[layer][iView]->SetPoint(iDiv, pos, peak);
	m_chargeStrip[layer][iView]->SetPointError(iDiv, errPos, errPeak);

	if( peak > 0.0 ){
	  float chargeScale = 5.0 / peak;
	  if( fabs(chargeScale-1) > 0.3 ){
	    std::cout << "WARNIN, Abnormal charge scale: " << chargeScale 
		      << ", (L,V,FE)=(" << layer << ", " << iView << ", " 
		      << iDiv << ")" << std::endl;
	    m_log << "WARNIN, Abnormal charge scale: " << chargeScale 
		  << ", (L,V,FE)=(" << layer << ", " << iView << ", " 
		  << iDiv << ")" << std::endl;
	    if( chargeScale > 1.3 ) chargeScale = 1.3;
	    if( chargeScale < 0.7 ) chargeScale = 0.7;
	  }
	  m_chargeScale[layer][iView][iDiv] = chargeScale;
	}
	
	m_log << "Charge " << layer << ' ' << iView << ' ' << pos << ' '
	       << area << ' ' << ave << ' ' << rms << ", " << *(par+0) << ' '
	       << *(par+1) << ' ' << *(par+2) << ' ' << *(par+3)
	       << std::endl;
	
      }
      cout << endl;
    }
  }
}


bool totCalib::readTotConvFile(const char* dir, const char* runid)
{
  string filename;
  char fname[] = "/398000364/TkrTotGainNt_LayerY17_398000364.tnt";
  for(int layer = 0; layer != g_nPlane; ++layer) {
    for(int iView = 0; iView != g_nView; ++iView) {
      filename = dir;
      if( iView == 0 )
	sprintf(fname,"/%s/TkrTotGainNt_LayerX%d_%s.tnt", runid, layer, runid);
      else
	sprintf(fname,"/%s/TkrTotGainNt_LayerY%d_%s.tnt", runid, layer, runid);
      filename += fname;
      if( !readTotConv( layer, iView, filename.c_str() ) )
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


bool totCalib::readTotConvXmlFile(const char* dir, const char* runid)
{
  string filename;
  filename = dir;
  char fname[] = "/398000364/TkrTotGain_398000364.xml";
  sprintf(fname,"/%s/TkrTotGain_%s.xml", runid, runid);

  filename += fname;
  std::cout << "Open xml file: " << filename << std::endl;
  m_log << "TOT xml file: " << filename << std::endl;
  


  xml::XmlParser* parser = new xml::XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse(filename.c_str());
  }
  catch (xml::ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    delete parser;
    return false;
  }
  if (doc != 0) {  // successful
    std::cout << "Document successfully parsed" << std::endl;

    // look up generic attributes
    DOMElement* docElt = doc->getDocumentElement();
    DOMElement* attElt = xml::Dom::findFirstChildByName(docElt,"generic");
    //DOMElement* isElt = xml::Dom::findFirstChildByName(attElt,"inputSample");
    
    try {
      m_tot_runid  = xml::Dom::getAttribute(attElt, "runId");
      //m_timeStamp = xml::Dom::getAttribute(attElt, "timestamp");
    }
    catch (xml::DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    // look up tower attributes
    attElt = xml::Dom::findFirstChildByName(docElt,"tower");
    try {
      m_tower_col  = xml::Dom::getIntAttribute(attElt, "col");
      m_tower_row = xml::Dom::getIntAttribute(attElt, "row");
      m_tower_serial  = xml::Dom::getAttribute(attElt, "hwserial");
    }
    catch (xml::DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    std::cout << "tower row: " << m_tower_row << ", col: " << m_tower_col 
	      << ", serial: " << m_tower_serial 
	      << ", runid: " << m_tot_runid << ", timeStamp: " << m_timeStamp
	      << std::endl;

    XMLCh* xmlchElt = XMLString::transcode("uniplane");
    DOMNodeList* conList = doc->getElementsByTagName(xmlchElt);
    int len = conList->getLength();   
    if( len != g_nPlane*g_nView ){
      std::cout << "ERROR: # of layers in xml is invalid, " << len << std::endl;
      m_log << "ERROR: # of layers in xml is invalid, " << len << std::endl;
      return false;
    }

    for(int i=0; i<len; i++){//each layers loop
      DOMNode* childNode = conList->item(i);
      int tray = xml::Dom::getIntAttribute(childNode, "tray");
      std::string which = xml::Dom::getAttribute(childNode, "which");
      std::cout << "(tray,which)=(" << tray << ", " << which << ") ";
     
      //get first child element
      DOMElement* elder = xml::Dom::getFirstChildElement(childNode);
      DOMElement* younger;

      int view = (tray+1) % 2;
      int layer = tray;
      if( which == "bot" ) layer -= 1;
      if( layer >= g_nPlane || layer < 0 ){
	std::cout << "Invalid layer id: " << layer << std::endl;
	m_log << "Invalid layer id: " << layer << std::endl;
	continue;
      }
	
      int numStrip = 0;
      while( getParam( elder, layer, view ) ){
	younger = xml::Dom::getSiblingElement( elder );
	elder = younger;
	numStrip++;
      }
      if( numStrip != g_nStrip ){
	std::cout << "ERROR: # of strips in xml is invalid, " 
		  << numStrip << std::endl;
	m_log << "ERROR: # of strips in xml is invalid, " 
	      << numStrip << std::endl;
	return false;
      }
    }
  }
  else return false;
  return true;
}

bool totCalib::getParam(const DOMElement* totElement, int layer, int view){  
  int stripId;
  double quad,gain,offset;
  try{
    stripId = xml::Dom::getIntAttribute( totElement, "id" );
  } //if there isn't next strip,go to next layer or view
  catch(xml::DomException ex){
    cout << "finished (layer,view)=(" << layer << ", "<< view << ")" << endl;
    return false;
  }
  if( stripId < 0 || stripId >= g_nStrip ){
    std::cout << "ERROR: (L,V)=(" << layer << ", " << view 
	      << "), Invalid strip id: " << stripId << std::endl;
    m_log << "ERROR: (L,V)=(" << layer << ", " << view 
	  << "), Invalid strip id: " << stripId << std::endl;
    return true;
  }

  try{
    quad = xml::Dom::getDoubleAttribute(totElement,"quad");
  }
  catch(xml::DomException ex){
    cout << "ERROR, no attribute for quad: (L,V,S)=(" << layer << ", " 
	 << view << ", "  << stripId << ")" << endl;
    m_log << "ERROR, no attribute for quad: (L,V,S)=(" << layer << ", " 
	  << view << ", "  << stripId << ")" << endl;
  }
  try{
    gain = xml::Dom::getDoubleAttribute(totElement,"slope");
  }
  catch(xml::DomException ex){
    cout << "ERROR, no attribute for slope: (L,V,S)=(" << layer << ", " 
	 << view << ", "  << stripId << ")" << endl;
    m_log << "ERROR, no attribute for slope: (L,V,S)=(" << layer << ", " 
	  << view << ", "  << stripId << ")" << endl;
  }
  try{
    offset = xml::Dom::getDoubleAttribute(totElement,"intercept");
  }
  catch(xml::DomException ex){
    cout << "ERROR, no attribute for offset: (L,V,S)=(" << layer << ", " 
	 << view << ", "  << stripId << ")" << endl;
    m_log << "ERROR, no attribute for offset: (L,V,S)=(" << layer << ", " 
	  << view << ", "  << stripId << ")" << endl;
  }
  /*  cout <<"stripId" << stripId
       <<",offset" << offset
       <<",gain" << gain
       <<",quad" << quad <<endl;*/
  m_totOffset[layer][view][stripId] = offset;
  m_totGain[layer][view][stripId] = gain;
  m_totQuadra[layer][view][stripId] = quad;
  return true;
}



float totCalib::calcCharge(int layer, int view, int iStrip, int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  // TOT to charge conversion
  float charge = m_totOffset[layer][view][iStrip] 
    + time*m_totGain[layer][view][iStrip]
    + time*time*m_totQuadra[layer][view][iStrip];
  
  return charge;
}


void totCalib::fillXml()//takuya
{

  std::string filename = m_outputDir;
  char fname[] = "/TkrFMX_TkrChargeScale_050131-161012.xml";
  
  sprintf( fname, "/%s_TkrChargeScale_%s-%s.xml", 
	   m_tower_serial.c_str(), m_dateStamp.c_str(), m_timeStamp.c_str() );

  filename += fname;

  std::ofstream output( filename.c_str() );
  if( output ){
    std::cout << "Open output xml file: " << filename << std::endl;
    m_log << "Output xml file: " << filename << std::endl;
  }
  else{
    std::cout << filename << " cannot be opened." << std::endl;
    return;
  }


  output << "<?xml version='1.0' ?>" << std::endl
	 << "<!DOCTYPE chargeScale [" << std::endl
	 << "<!-- dtd for Tracker calibrations -->" << endl
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
	 <<"<!ATTLIST inputSample  startTime CDATA    #REQUIRED"<<endl
	 <<"                       stopTime  CDATA    #REQUIRED"<<endl
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
	 <<"   data_thr  NMTOKEN #IMPLIED >"<<endl
	 << "]>" << std::endl;


  output << "<chargeScale>" << endl
	 << "   <generic calType='ChargeScale' creatorName='totCalib'"
	 << " creatorVersion ='" << m_version
	 << "' fmtVersion='NA' instrument='TWR' runId='" << m_tot_runid 
	 << "' timestamp='" << m_dateStamp << m_timeStamp << "'>" << std::endl
	 << "    <inputSample mode='NA' source='CosmicMuon' startTime='" 
	 << m_startTime << "' stopTime='" << m_stopTime 
	 << "' triggers='TKR'>" << std::endl
	 << " Cosmic ray muon data for charge scale calibration " << std::endl
	 << "    </inputSample>" << std::endl
	 << "  </generic>" << std::endl
	 << "  <tower row='" << m_tower_row << "' col='" << m_tower_col 
	 << "' hwserial='" << m_tower_serial << "'>" << endl;

  for(int layer = 0; layer != g_nPlane; ++layer) {
    for(int iView = 0; iView != g_nView; ++iView) {
      int tray;
      string which;
      if(iView==0){
	tray = 2 * ( layer/2 ) + 1;
	if(layer%2==0) which = "bot";
	else which = "top";
      }
      else{
	tray = 2 * ( (layer+1)/2 );
	if(layer%2==0) which = "top";
	else which = "bot";
      }
      if( iView == 0 )
	output << "   <!-- layer X" << layer << " -->" << std::endl;
      else
	output << "   <!-- layer Y" << layer << " -->" << std::endl;
      output << "   <uniplane tray='" << tray << "' which='"
	     << which << "'>" << endl;
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	output << "    <gtfeScale id='" << iDiv << "' chargeScale='" 
               <<  m_chargeScale[layer][iView][iDiv] << "'/>" << endl;
      }
      output << "   </uniplane>" << endl; 
    }
  }
  output     << "  </tower>" << endl;  
  output     << " </chargeScale>" << endl;
}


//-----------------------------------------------------------------------
//
//      Convoluted Landau and Gaussian Fitting Funcion
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
  if( sigma<=0.0 || width<=0.0 ) return 0.0;

  // MP shift correction
  mpc = par[1] - mpshift * width;
  if( mpc < 0.0 ) return 0.0;
  
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
