
#include <cmath>
#include <ctime>
#include <cassert>

#include "totCalib.h"
#include "facilities/Util.h"
 
using std::string;
using std::cout;
using std::endl;

XERCES_CPP_NAMESPACE_USE


totCalib::totCalib( const std::string analysisType = "MIP calibration" ): 
  m_reconFile(0), m_reconTree(0), 
  m_reconEvent(0), m_digiFile(0), m_digiTree(0),
  m_digiEvent(0), m_totFile(0)
{
  
  std::string m_dtd = "$(CALIBUTILROOT)/xml/";
  int status = facilities::Util::expandEnvVar(&m_dtd);
  if(status==-1)
    {
      std::cout<<m_dtd.c_str()<<" not found!"<<std::endl;
    }

  std::cout << analysisType << std::endl;
  if( analysisType == "badStrips" ) 
    {
      m_badStrips = true;
      m_dtd += "badStrips.dtd";
    }
  else {
    m_badStrips = false;
    m_dtd += "tkrCalibration.dtd";
  }
  std::cout << "DTD file: " << m_dtd << std::endl;

  // get version number from CVS string
  std::string tag = "$Name:  $";
  int i = tag.find( " " );
  tag.assign( tag, i+1, tag.size() );
  i = tag.find( " " );
  tag.assign( tag, 0, i ) ;
  m_tag = tag;

  std::string version = "$Revision: 1.23 $";
  i = version.find( " " );
  version.assign( version, i+1, version.size() );
  i = version.find( " " );
  version.assign( version, 0, i ) ;
  m_version = version;
  std::cout << "Tag: " << m_tag << ", version: " << m_version << std::endl;

  m_first_run = 999999999;
  m_last_run = 0;
  m_startTime="01/20 2005, 22:52 GMT";
  m_stopTime="01/20 2005, 22:52 GMT";

  // get time 
  time_t rawtime=0;
  time( &rawtime );
  // format time into a null-terminated string
  char nts[] = "050124-000000";
  size_t ntsmax=10;
  strftime( nts, ntsmax, "%y%m%d", gmtime( &rawtime ) );
  m_dateStamp = nts;
  strftime( nts, ntsmax, "%H%M%S", gmtime( &rawtime ) );
  m_timeStamp = nts;

  if( m_badStrips ){
    m_occDist = new TH1F("occDist", "occDist", 200, 0, 200);
    m_poissonDist = new TH1F("poissonDist", "poissonDist", 40, -20, 0);
    m_aPos[0] = new TH1F("apos0", "apos0", 100, -50, 50);
    m_aPos[1] = new TH1F("apos1", "apos1", 100, -50, 50);
    m_aPos[2] = new TH1F("apos2", "apos2", 100, -50, 50);
    m_aPos[3] = new TH1F("apos3", "apos3", 100, -50, 50);
  }

  for(int tower = 0; tower != g_nTower; ++tower) {
    m_tower_serial[tower] = "None";
    m_tot_runid[tower] = "-1";

    for(int layer = 0; layer != g_nLayer; ++layer) {
      for(int iView = 0; iView != g_nView; ++iView) {
	char vw = 'X';
	if( iView != 0 ) vw = 'Y';
	if( !m_badStrips ){
	  char name[] = "var00000";
	  sprintf(name,"var%d%d%d", tower, layer, iView);
	  m_totStrip[tower][layer][iView] = new TGraphErrors(g_nDiv);
	  m_totStrip[tower][layer][iView]->SetName(name);
	  
	  char temp[] = "varCorr00000";
	  sprintf(temp,"varCorr%d%d%d", tower, layer, iView);
	  m_chargeStrip[tower][layer][iView] = new TGraphErrors(g_nDiv);
	  m_chargeStrip[tower][layer][iView]->SetName(temp);
	}
	if( m_badStrips ){
	  for(int iWafer = 0; iWafer != g_nWafer; ++iWafer) {
	    char name1[] = "occT00X17w3";
	    sprintf(name1,"occT%d%c%dw%d", tower, vw, layer, iWafer);
	    m_nHits[tower][layer][iView][iWafer] = new TH1F(name1, name1, 1536, 0, 1536);
	  }
	}
	else{
	  for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	    char name1[] = "totT00X17fe0004";
	    sprintf(name1,"totT%d%c%dfe%d", tower, vw, layer, iDiv);
	    m_totHist[tower][layer][iView][iDiv] = new TH1F(name1, name1, 100, 0, 200);
	    
	    char name2[] = "chargeT00X00fe0000";
	    sprintf(name2,"chargeT%d%c%dfe%d", tower, vw, layer, iDiv);
	    m_chargeHist[tower][layer][iView][iDiv] = new TH1F(name2, name2, 200, 0, 20);
	  }
	}
      }
    }
  }
}

totCalib::~totCalib() 
{
  if(m_totFile == 0) return;

  std::cout << "save histgrams" << std::endl;

  m_totFile->cd();

  if( m_badStrips ){
    m_occDist->Write(0, TObject::kOverwrite);
    m_poissonDist->Write(0, TObject::kOverwrite);
    m_aPos[0]->Write(0, TObject::kOverwrite);
    m_aPos[1]->Write(0, TObject::kOverwrite);
    m_aPos[2]->Write(0, TObject::kOverwrite);
    m_aPos[3]->Write(0, TObject::kOverwrite);
  }

  for(int tower = 0; tower != g_nTower; ++tower) {
    for(int layer = 0; layer != g_nLayer; ++layer) {
      
      for(int iView = 0; iView != g_nView; ++iView) {
	
	if( m_badStrips )
	  for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
	    m_nHits[tower][layer][iView][iWafer]->Write(0, TObject::kOverwrite);
	else{
	  m_totStrip[tower][layer][iView]->Write(0, TObject::kOverwrite);
	  m_chargeStrip[tower][layer][iView]->Write(0, TObject::kOverwrite);
	  for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	    m_totHist[tower][layer][iView][iDiv]->Write(0, TObject::kOverwrite);
	    m_chargeHist[tower][layer][iView][iDiv]->Write(0, TObject::kOverwrite);
	  }
	}
	
      }
    }
  }

  m_totFile->Close();
}

bool totCalib::setOutputFiles( const char* outputDir )
{
  m_outputDir = outputDir;
  std::string testId = "TE603";
  if( m_badStrips ) testId = "TE403";

  std::string filename;
  char fname[] = "/TE603_050121-000000.root";

  filename = m_outputDir;
  sprintf( fname, "/%s_%s-%s.log", testId.c_str(), 
	   m_dateStamp.c_str(), m_timeStamp.c_str() );
  filename += fname;
  m_log.open( filename.c_str() );
  if( m_log )
    std::cout << "Open log file: " << filename << std::endl;
  else{
    std::cout << filename << " cannot be opened." << std::endl;
    return false;
  }

  filename = m_outputDir;
  sprintf( fname, "/%s_%s-%s.root", testId.c_str(), 
	   m_dateStamp.c_str(), m_timeStamp.c_str() );
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


int totCalib::setInputRootFiles( const char* rootDir, const char* digiPrefix,
				 const char* reconPrefix, 
				 const std::vector<std::string>& runIds )
{

  std::string digiFile, reconFile;
  char fname[] = "135000933/v4r060302p8/calib-v1r0/grRoot/recon-v3r1p2_135000933_recon_RECON_100.root";
  TChain* digiChain = new TChain("Digi");
  TChain* reconChain = new TChain("Recon");

  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {

    int runid = atoi( (*run).c_str() );
    if( runid < m_first_run ) m_first_run = runid;
    if( runid > m_last_run ) m_last_run = runid;

    digiFile = rootDir;
    sprintf(fname,"/%d/%s_%d_digi_DIGI.root",
	    runid, digiPrefix, runid);
    digiFile += fname;
    std::cout << "open digi file: " << digiFile << endl;
    m_log << "digi file: " << digiFile << endl;
    digiChain->Add( digiFile.c_str() );

    int split = 0;
    while( true ){
      reconFile = rootDir;
      if( split == 0 )
	sprintf(fname,"/%d/%s_%d_recon_RECON.root",
		runid, reconPrefix, runid);
      else
	sprintf(fname,"/%d/%s_%d_recon_RECON_%d.root",
		runid, reconPrefix, runid, split);
      reconFile += fname;
      if( ! checkFile( reconFile ) ){
	if( split == 0 ) return false;
	else break;
      }
      std::cout << "open recon file: " << reconFile << endl;
      m_log << "recon file: " << reconFile << endl;
      reconChain->Add( reconFile.c_str() );
      split++;
    }

  }
  return setInputRootFiles( digiChain, reconChain );

}

bool totCalib::checkFile( const std::string fname ){

  bool flag;
  ifstream file( fname.c_str(), std::ios::in | std::ios::binary );
  flag = file.good();
  file.clear();
  file.close();

  return flag;
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
#ifdef OLD_RECON
  using namespace xml;
#else
  using namespace xmlBase;
#endif

  XmlParser* parsercReport = new XmlParser(true);
  DOMDocument* docrcReport = 0;
  try{
    docrcReport = parsercReport -> parse(reportFile);
  }
  catch (ParseException ex) {
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
      timeStamp = Dom::getAttribute(rcElt, "timestamp");
    }
    catch (DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    std::vector<std::string> keywords, values;
    keywords.push_back("RunId");
    keywords.push_back("StartTime");
    keywords.push_back("EndTime");
    keywords.push_back("SerialNos");

    for( unsigned int i=0; i<keywords.size(); i++){
      DOMElement* childElt 
	= Dom::findFirstChildByName( rcElt, keywords[i].c_str() );
      try {
	values.push_back( Dom::getTextContent(childElt) );
	//std::cout << keywords[i] << ": " << values[i] << std::endl;
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
    }

    std::string serials = values[3], towerId, tower_serial;
    while( true ){
      int pos = serials.find( "GTEM" );
      if( pos == string::npos ) break;
      serials.assign( serials, pos+5, serials.size() );
      pos = serials.find( "," );
      towerId.assign( serials, 0, pos );
      pos = serials.find( "tkr" );
      serials.assign( serials, pos+5, serials.size() );
      pos = serials.find( "'" );
      serials.assign( serials, pos+1, serials.size() );
      pos = serials.find( "'" );
      tower_serial.assign( serials, 0, pos );
      if( pos == string::npos ) break;
      int tower = atoi( towerId.c_str() );
      if( tower >= 0 && tower < g_nTower ){
	if( m_tower_serial[ tower ] == "None" ){
	  std::cout << towerId << " " << tower_serial << std::endl;
	  m_log << towerId << " " << tower_serial << std::endl;
	  m_towerList.push_back( tower );
	  m_tower_serial[ tower ] = tower_serial;
	}
	else if( m_tower_serial[ tower ] != tower_serial ){
	  std::cout << "Inconsistent tower serial IDs for tower " << tower
		    << ": " << m_tower_serial[ tower ] << " " << tower_serial
		    << std::endl;
	  m_log << "Inconsistent tower serial IDs for tower " << tower
		<< ": " << m_tower_serial[ tower ] << " " << tower_serial
		<< std::endl;
	  return false;
	}
      }
      else {
	std::cout << "Invalid tower number, contents of SerialNos:" 
		  << values[3] << std::endl;
	m_log << "Invalid tower number, contents of SerialNos:" 
	      << values[3] << std::endl;
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

  if( m_badStrips ){
    findBadStrips( nEvents );
    fillBadStrips();
  }
  else{
    fitTot();
    fillXml();//takuya
  }
}

void totCalib::analyzeEvent(int nEvents) 
{
  int mEvent = nEvents * 0.01;
  if( mEvent < 100 ) mEvent = 100;
  time_t startTime, currentTime;
  time( &startTime );

  for(int iEvent = 0; iEvent != nEvents; ++iEvent) {    
    if( iEvent >= mEvent ){
      time( &currentTime );
      int elapsedTime = currentTime - startTime;
      if( elapsedTime <= 0 ) elapsedTime = 1;
      int rEvents = nEvents - iEvent;
      std::cout << "# of events: " << iEvent << " (" << iEvent*101/nEvents 
		<< "%) in " << elapsedTime << " s, "
		<< iEvent/elapsedTime << " events/s, "
		<< rEvents << " events, "
		<< int(1.0*rEvents*elapsedTime/iEvent)
		<< " s to go" << std::endl;
      if( mEvent > nEvents*0.095 ) mEvent += nEvents * 0.1;
      else mEvent += nEvents * 0.01;
    }
    if( m_reconEvent ) m_reconEvent->Clear();
    if( m_digiEvent ) m_digiEvent->Clear();

    m_reconTree->GetEntry(iEvent);
    m_digiTree->GetEntry(iEvent);

    assert(m_reconEvent != 0);
    assert(m_digiEvent != 0);

    if(! passCut()) continue;

    if( m_badStrips ) fillOccupancy();
    else{
      getTot();
      fillTot();
    }
  }
  std::cout << "Data scan finished." << std::endl;
  time( &currentTime );
  std::cout << "total # of events: " << nEvents 
		<< " in " << (currentTime-startTime) << " s, "
		<< nEvents/(currentTime-startTime) << " events/s"
		<< std::endl;
  m_log << "total # of events: " << nEvents 
	<< " in " << (currentTime-startTime) << " s, "
	<< nEvents/(currentTime-startTime) << " events/s"
	<< std::endl;

}

void totCalib::getTot()
{
  int noOfTkrDigis = m_digiEvent->getTkrDigiCol()->GetLast()+1;

  for(int i = 0; i != noOfTkrDigis; ++i) {
    const TkrDigi* tkrDigi = m_digiEvent->getTkrDigi(i);

    assert(tkrDigi != 0);
    int iTower = tkrDigi->getTower().id();
    int iLayer = tkrDigi->getBilayer();

    GlastAxis::axis viewId = tkrDigi->getView();
    int view = (viewId == GlastAxis::X) ? 0 : 1;
    m_tot[iTower][iLayer][view][0] = tkrDigi->getToT(0);
    m_tot[iTower][iLayer][view][1] = tkrDigi->getToT(1);
    m_lastRC0Strip[iTower][iLayer][view] = tkrDigi->getLastController0Strip();
  } 
}
 
void totCalib::retrieveCluster()
{
#ifdef OLD_RECON
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
 
  assert(tkrRecon != 0);

  TObjArray* siClusterCol = tkrRecon->getClusterCol();

  int noOfTkrClusters = siClusterCol->GetLast()+1;
  for(int i = 0; i != noOfTkrClusters; ++i) {
    TkrCluster* cluster = dynamic_cast<TkrCluster*>(siClusterCol->At(i));
    m_cluster[cluster->getId()] = cluster;
  }
#endif
}

int totCalib::findTot(int towerId, int layerId, int view , int stripId)
{
  if(stripId <= m_lastRC0Strip[towerId][layerId][view] )
    return m_tot[towerId][layerId][view][0];
  else
    return m_tot[towerId][layerId][view][1];
}

void totCalib::fillTot() 
{
  retrieveCluster();

  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();

#ifdef OLD_RECON
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

#ifdef OLD_RECON
    int tower = cluster->getTower();
#else
    int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
#endif

    int layer = g_nLayer - planeId - 1;
    int view = (viewId == TkrCluster::X) ? 0 : 1;

    // require only a single strip
    if(cluster->getSize() != 1) continue;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip) {

      int tot = findTot(tower,layer, view, iStrip);
#ifndef OLD_RECON
      tot = cluster->getRawToT(); // we might want to check that this one is correct
#endif
      if( tot == 0 ) continue;

      float charge = calcCharge(tower,layer, view, iStrip, tot);

      static int nStripPerGroup = g_nStrip / g_nDiv;

      m_totHist[tower][layer][view][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
      m_chargeHist[tower][layer][view][iStrip/nStripPerGroup]->Fill(charge*(-m_dir.z()));
    }
  }
#else
  TObjArray* tracks = tkrRecon->getTrackCol();
  TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->First());

  TIter trk1HitsItr(tkrTrack);
  TkrTrackHit* pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {
    
    const TkrCluster* cluster = pTrk1Hit->getClusterPtr();
    if(cluster) {
      //a cluster is attached to the hit: proceed
      
      int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY()).id();
      int layer = cluster->getLayer();
      int view = cluster->getTkrId().getView();

      // require only a single strip
      if(cluster->getSize() != 1) continue;

      for(int iStrip = cluster->getFirstStrip(); 
	  iStrip != int(cluster->getLastStrip()+1); ++iStrip) {
	
	int tot = cluster->getRawToT();
	if( tot == 0 ) continue;
	
	//	std::cout<<"plane: "<<layer<<" view: "<<view<<" digi tot: "<<tot<<std::endl;
	//	std::cout<<m_totX[planeId][0]<<" "<<m_totX[planeId][1]<<" "<<m_totY[planeId][0]<<" "<<m_totY[planeId][1]<<std::endl;
	
	float charge = calcCharge(tower,layer, view, iStrip, tot);
	
	static int nStripPerGroup = g_nStrip / g_nDiv;
	
	m_totHist[tower][layer][view][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
	m_chargeHist[tower][layer][view][iStrip/nStripPerGroup]->Fill(charge*(-m_dir.z()));
      }
    }
  }
#endif
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
  
  char cvw[] = "XY";
  for(int tower = 0; tower != g_nLayer; ++tower){
    for(int layer = 0; layer != g_nLayer; ++layer) {
      for(int iView = 0; iView != g_nView; ++iView) {
	for(int iDiv = 0; iDiv != g_nDiv; ++iDiv){
	  std::cout << "Layer: " << cvw[iView] << layer 
		    << ", FE: " << iDiv << std::endl;
	  m_chargeScale[tower][layer][iView][iDiv] = 1.0;
	  
	  // fit uncorrected tot for each strip
	  float area = m_totHist[tower][layer][iView][iDiv]->GetEntries();
	  float ave = m_totHist[tower][layer][iView][iDiv]->GetMean();
	  float rms = m_totHist[tower][layer][iView][iDiv]->GetRMS();
	  //std::cout << area << " " << ave << " " << rms << std::endl;
	  if( area<100 || ave==0.0 || rms==0.0 ){ 
	    std::cout << "Layer: " << cvw[iView] << layer
		      << ", FE: " << iDiv << ", Entries: " << area
		      << ", Mean: " << ave << ", RMS: " << rms 
		      << " skipped." << std::endl;
	    m_log << "Layer: " << cvw[iView] << layer
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
	  
	  m_totStrip[tower][layer][iView]->SetPoint(iDiv, pos, peak);
	  m_totStrip[tower][layer][iView]->SetPointError(iDiv, errPos, errPeak);
	  
	  /*m_log << "Uncorrected tot " << layer << ' ' << iView << ' ' << pos
	    << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	    << errWidth << endl; */
	  
	  // fit charge for each strip
	  area = m_chargeHist[tower][layer][iView][iDiv]->GetEntries();
	  ave = m_chargeHist[tower][layer][iView][iDiv]->GetMean();
	  rms = m_chargeHist[tower][layer][iView][iDiv]->GetRMS();
	  //std::cout << area << " " << ave << " " << rms << std::endl;
	  if( area<100 || ave==0.0 || rms==0.0 ){ 
	    m_log << "Layer: " << cvw[iView] << layer
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
	  m_chargeHist[tower][layer][iView][iDiv]->Fit( "langau", "RBQ" );
	  
	  //0:width(scale) 1:peak 2:total area 3:width(sigma)
	  par = ffit->GetParameters();
	  error = ffit->GetParErrors();
	  //par = (m_chargeHist[layer][iView][iDiv]->GetFunction("landau"))->GetParameters();
	  //error = (m_chargeHist[layer][iView][iDiv]->GetFunction("landau"))->GetParErrors();
	  
	  peak = float( *(par+1) );
	  errPeak = float( *(error+1) );
	  
	  width = float( *(par+3) );
	  errWidth = float( *(error+3) );
	  
	  m_chargeStrip[tower][layer][iView]->SetPoint(iDiv, pos, peak);
	  m_chargeStrip[tower][layer][iView]->SetPointError(iDiv, errPos, errPeak);
	  
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
	    m_chargeScale[tower][layer][iView][iDiv] = chargeScale;
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
  

#ifdef OLD_RECON
  //typedef xml xmlBase;
  using namespace xml;
#else
  using namespace xmlBase;
#endif

  XmlParser* parser = new XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse(filename.c_str());
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    delete parser;
    return false;
  }
  if (doc != 0) {  // successful
    std::cout << "Document successfully parsed" << std::endl;

    // look up generic attributes
    DOMElement* docElt = doc->getDocumentElement();
    DOMElement* attElt = Dom::findFirstChildByName(docElt,"generic");
    //DOMElement* isElt = xml::Dom::findFirstChildByName(attElt,"inputSample");
    
    std::string tot_runid, hwserial;
    try {
      tot_runid  = Dom::getAttribute(attElt, "runId");
      //m_timeStamp = xml::Dom::getAttribute(attElt, "timestamp");
    }
    catch (DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    // look up tower attributes
    attElt = Dom::findFirstChildByName(docElt,"tower");
    try {
      hwserial  = Dom::getAttribute(attElt, "hwserial");
    }
    catch (DomException ex) {
      std::cout << "DomException:  " << ex.getMsg() << std::endl;
    }

    int towerId=-1;
    for( unsigned int tw=0; tw<m_towerList.size(); tw++)
      if( m_tower_serial[ m_towerList[tw] ] == hwserial )
	towerId = m_towerList[ tw ];
    std::cout << "tower: " << towerId << ", serial: " << hwserial 
	      << ", runid: " << tot_runid << ", timeStamp: " << m_timeStamp
	      << std::endl;
    if( towerId < 0 ) return false;

    XMLCh* xmlchElt = XMLString::transcode("uniplane");
    DOMNodeList* conList = doc->getElementsByTagName(xmlchElt);
    int len = conList->getLength();   
    if( len != g_nLayer*g_nView ){
      std::cout << "ERROR: # of layers in xml is invalid, " << len << std::endl;
      m_log << "ERROR: # of layers in xml is invalid, " << len << std::endl;
      return false;
    }

    for(int i=0; i<len; i++){//each layers loop
      DOMNode* childNode = conList->item(i);
      int tray = Dom::getIntAttribute(childNode, "tray");
      std::string which = Dom::getAttribute(childNode, "which");
      std::cout << "(tray,which)=(" << tray << ", " << which << ") ";
     
      //get first child element
      DOMElement* elder = Dom::getFirstChildElement(childNode);
      DOMElement* younger;

      int view = (tray+1) % 2;
      int layer = tray;
      if( which == "bot" ) layer -= 1;
      if( layer >= g_nLayer || layer < 0 ){
	std::cout << "Invalid layer id: " << layer << std::endl;
	m_log << "Invalid layer id: " << layer << std::endl;
	continue;
      }
	
      int numStrip = 0;
      while( getParam( elder, towerId, layer, view ) ){
	younger = Dom::getSiblingElement( elder );
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

bool totCalib::getParam(const DOMElement* totElement, int tower, int layer, int view){  
#ifdef OLD_RECON
  //typdef xml xmlBase;
  using namespace xml;
#else
  using namespace xmlBase;
#endif

  int stripId;
  std::vector<std::string> keywords;
  std::vector<float> values;
  keywords.push_back( "intercept" );
  keywords.push_back( "slope" );
  keywords.push_back( "quad" );

  try{
    stripId = Dom::getIntAttribute( totElement, "id" ); 
  } //if there isn't next strip,go to next layer or view
  catch(DomException ex){
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

  for( unsigned int i=0; i<keywords.size(); i++){
    try{
      float value = Dom::getDoubleAttribute( totElement, keywords[i].c_str() );
      values.push_back( value );
    }
    catch(DomException ex){
      cout << "ERROR, no attribute for " << keywords[i] << ": (L,V,S)=(" 
	   << layer << ", " << view << ", "  << stripId << ")" << endl;
      m_log << "ERROR, no attribute for " << keywords[i] << ": (L,V,S)=(" 
	    << layer << ", " << view << ", "  << stripId << ")" << endl;
    }
  }
  /*  cout <<"stripId" << stripId
      <<",offset" << offset
      <<",gain" << gain
      <<",quad" << quad <<endl;*/
  
  m_totOffset[tower][layer][view][stripId] = values[0];
  m_totGain[tower][layer][view][stripId] = values[1];
  m_totQuadra[tower][layer][view][stripId] = values[2];
  return true;
}


float totCalib::calcCharge(int tower, int layer, int view, int iStrip, int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  // TOT to charge conversion
  float charge = m_totOffset[tower][layer][view][iStrip] 
    + time*m_totGain[tower][layer][view][iStrip]
    + time*time*m_totQuadra[tower][layer][view][iStrip];
  
  return charge;
}


void totCalib::fillXml()//takuya
{

  std::string filename = m_outputDir;
  char fname[] = "/TkrFMX_TkrChargeScale_050131-161012.xml";
  
  std::string tot_runid = m_tot_runid[ m_towerList[0] ];
  if( m_towerList.size() == 1 ){
    sprintf( fname, "/%s_TkrChargeScale_%s-%s.xml", 
	     m_tower_serial[ m_towerList[0] ].c_str(), 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
  }
  else{
    sprintf( fname, "/LAT_TkrChargeScale_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    for( unsigned int i=1; i<m_towerList.size(); i++){
      tot_runid += ':' + m_tot_runid[ m_towerList[i] ];
    }
  }
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
  std::ifstream dtd( m_dtd.c_str() );
  if( dtd ){
    std::cout << "Open dtd file: " << m_dtd << std::endl;
    m_log << "dtd file: " << m_dtd << std::endl;
  }
  else{
    std::cout << m_dtd << " cannot be opened." << std::endl;
    return;
  }

  output << "<?xml version=\"1.0\" ?>" << std::endl
	 << "<!DOCTYPE chargeScale [" << std::endl;

  std::string line;
  while( dtd ){
    getline(dtd, line);
    output << line << std::endl;
  }

  output << "]>" << std::endl;


  output << "<chargeScale>" << endl
	 << "   <generic calType=\"ChargeScale\" creatorName=\"totCalib\""
	 << " creatorVersion =\"" << m_version
	 << "\" fmtVersion=\"NA\" instrument=\"TWR\" runId=\"" << tot_runid 
	 << "\" timestamp=\"" << m_dateStamp << m_timeStamp << "\">" << std::endl
	 << "    <inputSample mode=\"NA\" source=\"CosmicMuon\" startTime=\"" 
	 << m_startTime << "\" stopTime=\"" << m_stopTime 
	 << "\" triggers=\"TKR\">" << std::endl
	 << " Cosmic ray muon data for charge scale calibration " << std::endl
	 << "    </inputSample>" << std::endl
	 << "  </generic>" << std::endl;

  output.precision(3);
  char cvw[] = "XY";

  for( unsigned int tw=0; tw<m_towerList.size(); tw++ ){
    int tower = m_towerList[ tw ];
    idents::TowerId twrId(tower); 
    int tower_row = twrId.ix();
    int tower_col = twrId.iy();
    std::string hwserial = m_tower_serial[ tower ];
    output << "  <tower row=\"" << tower_row << "\" col=\"" << tower_col 
	   << "\" hwserial=\"" << hwserial << "\">" << endl;

    for(int layer = 0; layer != g_nLayer; ++layer) {
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
	output << std::endl
	       << "   <!-- **** layer " << cvw[iView]  << layer 
	       << " **** -->" << std::endl
	       << "   <uniplane tray=\"" << tray << "\" which=\""
	       << which << "\">" << std::endl;
	for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	  output << "    <gtfeScale id=\"" << iDiv << "\" chargeScale=\"" 
		 <<  m_chargeScale[tower][layer][iView][iDiv] << "\"/>" << endl;
	}
	output << "   </uniplane>" << endl; 
      }
    }
    output     << "  </tower>" << endl;  
  }
  output     << " </chargeScale>" << endl;

}
  

void totCalib::fillOccupancy() 
{

  //initialize container
  int nHits[g_nTower][g_nLayer][g_nView][g_nWafer+1];
  for( unsigned int tw=0; tw<m_towerList.size(); tw++){
    int tower = m_towerList[tw];
    for( int layer=0; layer<g_nLayer; layer++)
      for( int view=0; view<g_nView; view++)
	for( int i=0; i<g_nWafer+1; i++) nHits[tower][layer][view][i] = 0;
  }

  retrieveCluster();

  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();

  TObjArray* tracks = tkrRecon->getTrackCol();
#ifdef OLD_RECON
  TkrKalFitTrack* tkrTrack = dynamic_cast<TkrKalFitTrack*>(tracks->At(0));

  int nHitPlane = tkrTrack->getNumHits();


  for(int iPlane = 0; iPlane != nHitPlane; ++iPlane) {
    const TkrHitPlane* plane = tkrTrack->getHitPlane(iPlane);
    std::map<int, TkrCluster*>::const_iterator itr = m_cluster.find(plane->getIdHit());
    assert(itr != m_cluster.end());
    TkrCluster* cluster = itr->second;
    int planeId = cluster->getPlane();
    TkrCluster::view viewId = cluster->getView();

#ifdef OLD_RECON
    int tower = cluster->getTower();
#else
    int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
#endif
    int layer = g_nLayer - planeId - 1;
    int view = (viewId == TkrCluster::X) ? 0 : 1;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip){
      nHits[tower][layer][view][iStrip/384]++;
      nHits[tower][layer][view][g_nWafer]++;
    }
  }

  bool display = false;
  float dx=0.0, dy=0.0, pos, apos;
  int view, aview;

  for(int iPlane = 0; iPlane != nHitPlane; ++iPlane) {

    const TkrHitPlane* plane = tkrTrack->getHitPlane(iPlane);
    std::map<int, TkrCluster*>::const_iterator itr = m_cluster.find(plane->getIdHit());
    assert(itr != m_cluster.end());
    TkrCluster* cluster = itr->second;
    int planeId = cluster->getPlane();
#ifdef OLD_RECON
    int tower = cluster->getTower();
#else
    int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
#endif
    TkrCluster::view viewId = cluster->getView();
    TVector3 position = cluster->getPosition();
    float deltax = m_pos.X()+m_dir.X()/m_dir.Z()*(position.Z()-m_pos.Z()) - position.X();
    float deltay = m_pos.Y()+m_dir.Y()/m_dir.Z()*(position.Z()-m_pos.Z()) - position.Y();

    int layer = g_nLayer - planeId - 1;
    if( viewId == TkrCluster::X ){
      view = 0;
      aview = 1;
      if( fabs( deltax - dx ) > 2.0  ){
	//std::cout << layer << " " << view << ", " << deltax << " " << dx
	//  << " ***********************" << std::endl;
	break;
      }
      deltax -= dx;
      dx += deltax;
      pos = deltax;
      apos = deltay;
    }
    else{
      view = 1;
      aview = 0;
      if( fabs( deltay - dy ) > 2.0  ){
	//std::cout << layer << " " << view << ", " << deltay << " " << dy
	//  << " ***********************" << std::endl;
	break;
      }
      deltay -= dy;
      dy += deltay;
      pos = deltay;
      apos = deltax;
    }

    //std::cout << layer << " " << view << ", " << pos << " " << apos
    //      << std::endl;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip)
      if( nHits[tower][layer][aview][g_nWafer] > 0 ){
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( nHits[tower][layer][aview][iw] > 0 ){
	    m_nHits[tower][layer][view][iw]->Fill( iStrip );
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
      }
      else
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( fabs( apos-89.5*(iw-1.5) ) < 42 ){
	    m_nHits[tower][layer][view][iw]->Fill( iStrip );
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
  }
  if( display ){
    for(int iPlane = 0; iPlane != nHitPlane; ++iPlane) {
      const TkrHitPlane* plane = tkrTrack->getHitPlane(iPlane);
      std::map<int, TkrCluster*>::const_iterator itr = m_cluster.find(plane->getIdHit());
      assert(itr != m_cluster.end());
      TkrCluster* cluster = itr->second;
      int planeId = cluster->getPlane();
#ifdef OLD_RECON
      int tower = cluster->getTower();
#else
      int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
#endif
      TkrCluster::view viewId = cluster->getView();
      
      int layer = g_nLayer - planeId - 1;
      int view = (viewId == TkrCluster::X) ? 0 : 1;

      std::cout << tower << " " << layer << " " << view;
      for(int iStrip = cluster->getFirstStrip(); 
	  iStrip != int(cluster->getLastStrip()+1); ++iStrip)
	std::cout << " " << iStrip;
      std::cout << std::endl;
    }
  }
#else
  TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->First());

//   int nHits[g_nLayer][g_nView][g_nWafer+1];
//   for( int layer=0; layer<g_nLayer; layer++)
//     for( int view=0; view<g_nView; view++)
//       for( int i=0; i<g_nWafer+1; i++) nHits[layer][view][i] = 0;

  TIter trk1HitsItr(tkrTrack);
  TkrTrackHit* pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {
    
    const TkrCluster* cluster = pTrk1Hit->getClusterPtr();
    if(cluster) 
      {
#ifdef OLD_RECON
	int tower = cluster->getTower();
	TkrCluster::view viewId = cluster->getView();
	int view = (viewId == TkrCluster::X) ? 0 : 1;
#else
	int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
	int  view = cluster->getTkrId().getView();
#endif
	int layer = cluster->getLayer();
	
	for(int iStrip = cluster->getFirstStrip(); 
	    iStrip != int(cluster->getLastStrip()+1); ++iStrip){
	  nHits[tower][layer][view][iStrip/384]++;
	  nHits[tower][layer][view][g_nWafer]++;
	}
    }
  }

  bool display = false;
  float dx=0.0, dy=0.0, pos, apos;
  int aview;

  trk1HitsItr.Reset();
  pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {
    
    const TkrCluster* cluster = pTrk1Hit->getClusterPtr();
    if(!cluster) continue;

#ifdef OLD_RECON
    int tower = cluster->getTower();
    TkrCluster::view viewId = cluster->getView();
    int view = (viewId == TkrCluster::X) ? 0 : 1;
#else
    int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
    int view = cluster->getTkrId().getView();
#endif
    int layer = cluster->getLayer();
    TVector3 position = cluster->getPosition();
    float deltax = m_pos.X()+m_dir.X()/m_dir.Z()*(position.Z()-m_pos.Z()) - position.X();
    float deltay = m_pos.Y()+m_dir.Y()/m_dir.Z()*(position.Z()-m_pos.Z()) - position.Y();

    if( view == 0 ){
      aview = 1;
      if( fabs( deltax - dx ) > 2.0  ){
	//std::cout << layer << " " << view << ", " << deltax << " " << dx
	//  << " ***********************" << std::endl;
	break;
      }
      deltax -= dx;
      dx += deltax;
      pos = deltax;
      apos = deltay;
    }
    else{
      aview = 0;
      if( fabs( deltay - dy ) > 2.0  ){
	//std::cout << layer << " " << view << ", " << deltay << " " << dy
	//  << " ***********************" << std::endl;
	break;
      }
      deltay -= dy;
      dy += deltay;
      pos = deltay;
      apos = deltax;
    }

    //std::cout << layer << " " << view << ", " << pos << " " << apos
    //      << std::endl;

    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip)
      if( nHits[layer][aview][g_nWafer] > 0 ){
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( nHits[tower][layer][aview][iw] > 0 ){
	    m_nHits[tower][layer][view][iw]->Fill( iStrip );
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
      }
      else
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( fabs( apos-89.5*(iw-1.5) ) < 42 ){
	    m_nHits[tower][layer][view][iw]->Fill( iStrip );
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
  }

  if( display ){
    trk1HitsItr.Reset();
    pTrk1Hit = 0;
    while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {      
      const TkrCluster* cluster = pTrk1Hit->getClusterPtr();
      if(!cluster) continue;

#ifdef OLD_RECON
      int tower = cluster->getTower();
      TkrCluster::view viewId = cluster->getView();
      int view = (viewId == TkrCluster::X) ? 0 : 1;
#else
      int tower = TowerId(cluster->getTkrId().getTowerX(),cluster->getTkrId().getTowerY())id();
      int view = cluster->getTkrId().getView();
#endif
      int layer = cluster->getLayer();
      
      //      int layer = g_nLayer - planeId - 1;
      //      int view = (viewId == TkrCluster::X) ? 0 : 1;

      std::cout << tower << " " << layer << " " << view;
      for(int iStrip = cluster->getFirstStrip(); 
	  iStrip != int(cluster->getLastStrip()+1); ++iStrip)
	std::cout << " " << iStrip;
      std::cout << std::endl;
    }
  }
#endif
}


void totCalib::findBadStrips( int nEvents )
{  

  float sumThreshold = sumThrPerEvent * nEvents;
  if( sumThreshold < 2.0 ) sumThreshold = 2.0;
  float occThreshold = occThrPerEvent * nEvents;
  if( occThreshold < 1.0 ) occThreshold = 1.0;

  double fac=1.0; //prepare factorial upto 99
  double factorial[200];
  factorial[0] = fac;
  for(int k=1;k!=200;k++){
    fac*=k;
    factorial[k] = fac;
  }

  for( unsigned int tw=0; tw!=m_towerList.size(); ++tw){
    int tower = m_towerList[ tw ]; 
    m_log << "Tower: " << tower << ", ID: " << m_tower_serial[tower] << std::endl;
    for(int layer = 0; layer != g_nLayer; ++layer) {
      for(int view = 0; view != g_nView; ++view) {
	char vw = 'X';
	if( view != 0 ) vw = 'Y';
	m_log << "Layer: " << vw << layer << std::endl;
	for( int strip=0; strip!=g_nStrip; strip++){
	  float sum = 0.0;
	  bool deadFlag = false;
	  int occ[g_nWafer];
	  for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	    int occupancy = m_nHits[tower][layer][view][iWafer]->GetBinContent( strip + 1 );
	    if(layer%2 == 0) occ[ g_nWafer-1-iWafer ] = occupancy;
	    else occ[ iWafer ] = occupancy;
	    sum += occupancy;
	    m_occDist->Fill( occupancy+0.1 );
	    if( occupancy < occThreshold ){
	      deadFlag = true;
	    }
	  }
	  
	  if( deadFlag || sum < sumThreshold ){
	    m_log << strip << ": " << sum << ", ";
	    int nBad = 0;
	    for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	      m_log << " " << occ[iWafer];
	    }
	    if( sum == 0 ){
	      //categorize as disconnected
	      m_deadStrips[tower][layer][view][1].push_back(strip); 
	      m_log << " *1*" << std::endl;
	      continue;
	    }
	    // intermittently disconnected
	    else if( sum < sumThreshold ) nBad = 3;
	    
	    sum = occ[0];
	    m_log << ", ";
	    bool flagBreak = false;
	    int occBreak = 0;
	    for(int iWafer = 1; iWafer != g_nWafer; ++iWafer){
	      int value = occ[iWafer];
	      sum += value;
	      double mean = sum / (iWafer+1);
	      if( value < 200 && mean < 199.5 ){
		float p_poisson = pow(mean,value)*exp(-mean)/factorial[value];
		float norm      = pow(mean,mean)*exp(-mean)/factorial[(int)(mean+0.5)];
		p_poisson = log( p_poisson/norm ); //normalize
		m_log << " " << (int)p_poisson;
		m_poissonDist->Fill( p_poisson );
		if( p_poisson < poissonThreshold ){ 
		  flagBreak = true;
		  occBreak += occ[iWafer];
		}
	      }
	    }
	    if( flagBreak )
	      if( occBreak == 0 ) nBad = 2; // partial disconnected
	      else nBad = 4; // intermittent partial disconnected
	    if( nBad > 0 ){
	      m_deadStrips[tower][layer][view][nBad].push_back(strip);
	      m_log << " *" << nBad << "*" << std::endl;
	      continue;
	    }
	    else m_log << std::endl;
	  }
	}
	
	m_log << vw << layer << ", # of bad channels:";
	std::cout << vw << layer << ", # of bad channel:";
	for( int iBad=0; iBad<g_nBad; iBad++){
	  m_log << " " << m_deadStrips[tower][layer][view][iBad].size();
	  std::cout << " " << m_deadStrips[tower][layer][view][iBad].size();
	}
	m_log << std::endl;
	std::cout << std::endl;
      }
    }
  }
}

void totCalib::fillBadStrips()
{

  std::string filename = m_outputDir;
  char fname[] = "/TkrFMX_DeadStrips_050131-161012.xml";
    
  if( m_towerList.size() == 1 )
    sprintf( fname, "/%s_DeadStrips_%s-%s.xml", 
	     m_tower_serial[ m_towerList[0] ].c_str(), 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
  else
    sprintf( fname, "/LAT_TkrChargeScale_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
  
  
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
  std::ifstream dtd( m_dtd.c_str() );
  if( dtd ){
    std::cout << "Open dtd file: " << m_dtd << std::endl;
    m_log << "dtd file: " << m_dtd << std::endl;
  }
  else{
    std::cout << m_dtd << " cannot be opened." << std::endl;
    return;
  }

  output << "<?xml version=\"1.0\" ?>" << std::endl
	 << "<!DOCTYPE badStrips [" << std::endl;
    
  std::string line;
  while( dtd ){
    getline(dtd, line);
    output << line << std::endl;
  }
    
  output << "]>" << std::endl;
    
    
  output << "<badStrips badType=\"dead\">" << std::endl
	 << "<!-- includes partial dead strips; "
	 << " intermitent and/or wrire bond broken " 
	 << "between SSD wafers -->" << std::endl;

  output << "  <generic calType=\"stripOccupancy\" creatorName=\"badStrips\""
	 << " creatorVersion =\"" << m_version
	 << "\" fmtVersion=\"NA\" instrument=\"TWR\" runId=\"" 
	 << m_first_run << '-' << m_last_run
	 << "\" timestamp=\"" << m_dateStamp << m_timeStamp << "\">" << std::endl
	 << "    <inputSample mode=\"NA\" source=\"CosmicMuon\" startTime=\"" 
	 << m_startTime << "\" stopTime=\"" << m_stopTime 
	 << "\" triggers=\"TKR\">" << std::endl
	 << " Cosmic ray muon data for occupancy analysis " << std::endl
	 << "    </inputSample>" << std::endl
	 << "  </generic>" << std::endl;

  //dead, disconnected, partial disconnected, intermittently disconnected, 
  //intermittently partial disconnexcted
  int howBad[g_nBad] = {2,4,12,20,28}; 
  std::string cBad[g_nBad] = {"dead","disconnected","partially disconnected",
			      "intermittently disconnected",
			      "intermittently partially connected"};
  char cvw[] = "XY";
  
  for( unsigned int tw=0; tw<m_towerList.size(); tw++ ){
    int tower = m_towerList[ tw ];
    idents::TowerId twrId(tower); 
    int tower_row = twrId.ix();
    int tower_col = twrId.iy();
    std::string hwserial = m_tower_serial[ tower ];
    output << "  <tower row=\"" << tower_row << "\" col=\"" << tower_col 
	   << "\" hwserial=\"" << hwserial << "\""
	   << " nOnbdCalib=\"false\" nOnbdTrig=\"false\""
	   << " nOnbdData=\"false\"" << ">" << std::endl;
    
    for(int layer = 0; layer != g_nLayer; ++layer) {
      for(int view = 0; view != g_nView; ++view) {
	int tray;
	string which;
	if(view==0){
	  tray = 2 * ( layer/2 ) + 1;
	  if(layer%2==0) which = "bot";
	  else which = "top";
	}
	else{
	  tray = 2 * ( (layer+1)/2 );
	  if(layer%2==0) which = "top";
	  else which = "bot";
	}
	
	output << std::endl
	       << "    <!-- layer " << cvw[view] << layer << " -->" << std::endl;
	
	for( int iBad=0; iBad!=g_nBad; iBad++ ){
	  int itr = m_deadStrips[tower][layer][view][iBad].size();
	  output << "    <!-- # of " << cBad[iBad] << " strips: " << itr 
		 << " -->" << std::endl 
		 << "    <uniplane tray=\"" << tray << "\" which=\""
		 << which << "\" nOnbdCalib=\"false\" nOnbdTrig=\"false\""
		 << " nOnbdData=\"false\" howBad=\"" << howBad[iBad] << "\"";
	  
	  if(itr){
	    output << ">" << std::endl << "      <stripList strips=\"";
	    for(int i=0;i!=itr;i++)
	      output << " " << m_deadStrips[tower][layer][view][iBad][i];
	    output << "\"/>" << std::endl 
		   << "    </uniplane>" << std::endl;
	  }
	  else output << "/>" << std::endl;
	  
	}
      }
    }
    output << "  </tower>" << std::endl;
  }
  output << "</badStrips>" << std::endl;
  output.close();
  dtd.close();
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
