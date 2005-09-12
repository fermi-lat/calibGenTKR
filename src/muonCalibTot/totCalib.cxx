
#include <cmath>
#include <ctime>
#include <cassert>

#include "totCalib.h"
#include "facilities/Util.h"
#ifndef OLD_RECON
#include "commonRootData/idents/TkrId.h"
#endif
#include "commonRootData/idents/TowerId.h"

using std::string;
using std::cout;
using std::endl;

XERCES_CPP_NAMESPACE_USE

//#define PRINT_DEBUG 1

//
// Cluster class implementation
//
Cluster::Cluster( int strip, int TOT=0 ){
  firstStrip = strip;
  lastStrip = strip;
  tower = -1;
  uniPlane = -1;
  tot = TOT;
  correctedTot = -100.0;
}

bool Cluster::addStrip( int strip ){

  if( strip == lastStrip+1 ){ // assume strip number is sorted.
    lastStrip = strip;
    return true;
  }
  else return false;

}

//
// layerId class implementation
//
inline layerId::layerId( int lyr, int vw, int twr ){ 
  setLayer( lyr, vw, twr ); }
inline layerId::layerId( int tr, std::string wh, int twr ){ 
  setTray( tr, wh, twr ); }
inline layerId::layerId( int unp ){ setUniPlane( unp ); }

inline void layerId::setLayer( int lyr, int vw, int twr ){
  layer = lyr; view = vw; tower = twr;
  layerToTray();
  trayToUniPlane();
}

inline void layerId::setUniPlane( int unp, int twr ){
  uniPlane = unp; tower = twr;
  uniPlaneToTray();
  trayToLayer();
}

inline void layerId::setTray( int tr, std::string wh, int twr ){
  tray = tr; which=wh; tower=twr;
  trayToUniPlane();
  trayToLayer();
}

inline void layerId::trayToUniPlane(){
  uniPlane = tray * 2;
  if( which == "bot" ) uniPlane--;
}

inline void layerId::trayToLayer(){
  view = (tray+1) % 2;
  layer = tray;
  if( which == "bot" ) layer--;
}

void layerId::layerToTray(){
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
}


inline void layerId::uniPlaneToTray(){
  tray = (uniPlane+1) / 2;
  if( uniPlane%2 == 0 ) which ="top";
  else which = "bot";
}

//
// towerVar class implementation
//
towerVar::towerVar( int twr, bool badStrips ){
  towerId = twr;
  hwserial = "None";
  runid = "-1";
  bsVar.clear();
  tcVar.clear();

  for( int unp=0; unp!=g_nUniPlane; unp++){
    badStripVar bsv;
    totCalibVar tcv;
    for( int strip=0; strip!=g_nStrip; strip++){
      rHits[unp][strip] = 0;
      dHits[unp][strip] = 0;
      if( badStrips ){
	bsv.tHits[strip] = 0;
	bsv.lHits[strip] = 0;
	for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
	  for( int tDiv = 0; tDiv != g_nTime; tDiv++)
	    bsv.nHits[strip][iWafer][tDiv] = 0;
      }
    }
    if( badStrips ) bsVar.push_back( bsv );
    else tcVar.push_back( tcv );
  }
#ifdef PRINT_DEBUG
  std::cout << "towerVar constructer " << twr << std::endl;
#endif
}

void towerVar::readHists( TFile* hfile, UInt_t iRoot, UInt_t nRoot ){

  TH1F* hist, *rhist, *dhist, *thist, *lhist;
  char name[] = "roccT00X17w3t4";
  char cvw[] = "XY";
  for( int unp=0; unp!=g_nUniPlane; unp++){
    layerId lid( unp );
    int layer = lid.layer;
    int view = lid.view;
    sprintf(name,"roccT%d%c%d", towerId, cvw[view], layer);
    rhist = (TH1F*)hfile->Get( name );
    sprintf(name,"doccT%d%c%d", towerId, cvw[view], layer);
    dhist = (TH1F*)hfile->Get( name );
    for( int strip=0; strip!=g_nStrip; strip++){      
      rHits[unp][strip] += (int)rhist->GetBinContent( strip+1 );
      dHits[unp][strip] += (int)dhist->GetBinContent( strip+1 );
    }
  }

  for( UInt_t unp=0; unp!=bsVar.size(); unp++){
    layerId lid( unp );
    int layer = lid.layer;
    int view = lid.view;
    sprintf(name,"toccT%d%c%d", towerId, cvw[view], layer);
    thist = (TH1F*)hfile->Get( name );
    sprintf(name,"loccT%d%c%d", towerId, cvw[view], layer);
    lhist = (TH1F*)hfile->Get( name );
    for( int strip=0; strip!=g_nStrip; strip++){
      bsVar[unp].lHits[strip] += (int)lhist->GetBinContent( strip+1 );
      bsVar[unp].tHits[strip] += (int)thist->GetBinContent( strip+1 );
    }
    for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
      for( int tDiv = 0; tDiv != g_nTime; tDiv++){
	sprintf(name,"occT%d%c%dw%dt%d", towerId, cvw[view], layer, iWafer, (tDiv+iRoot*g_nTime)/nRoot);
	hist = (TH1F*)hfile->Get( name );
	for( int strip=0; strip!=g_nStrip; strip++)
	  bsVar[unp].nHits[strip][iWafer][tDiv] 
	    += (int)hist->GetBinContent( strip+1 );
      }
  }
}


void towerVar::saveHists(){

  TH1F* hist, *rhist, *dhist, *thist, *lhist;
  char name[] = "roccT00X17w3t4";
  char cvw[] = "XY";
  for( int unp=0; unp!=g_nUniPlane; unp++){
    layerId lid( unp );
    int layer = lid.layer;
    int view = lid.view;
    sprintf(name,"roccT%d%c%d", towerId, cvw[view], layer);
    rhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    sprintf(name,"doccT%d%c%d", towerId, cvw[view], layer);
    dhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    for( int strip=0; strip!=g_nStrip; strip++){
      rhist->Fill( strip+0.1, rHits[unp][strip] );
      dhist->Fill( strip+0.1, dHits[unp][strip] );
    }
    rhist->Write(0, TObject::kOverwrite);
    dhist->Write(0, TObject::kOverwrite);
  }

  for( UInt_t unp=0; unp!=bsVar.size(); unp++){
    layerId lid( unp );
    int layer = lid.layer;
    int view = lid.view;
    sprintf(name,"toccT%d%c%d", towerId, cvw[view], layer);
    thist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    sprintf(name,"loccT%d%c%d", towerId, cvw[view], layer);
    lhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    for( int strip=0; strip!=g_nStrip; strip++){
      lhist->Fill( strip+0.1, bsVar[unp].lHits[strip] );
      thist->Fill( strip+0.1, bsVar[unp].tHits[strip] );
    }
    thist->Write(0, TObject::kOverwrite);
    lhist->Write(0, TObject::kOverwrite);
    for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
      for( int tDiv = 0; tDiv != g_nTime; tDiv++){
	sprintf(name,"occT%d%c%dw%dt%d", towerId, cvw[view], layer, iWafer, tDiv);
	hist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
	for( int strip=0; strip!=g_nStrip; strip++)
	  hist->Fill( strip+0.1, bsVar[unp].nHits[strip][iWafer][tDiv] );
	hist->Write(0, TObject::kOverwrite);
      }
  }
}

//
//
//
/*
histCol::histCol( int towerId ){
  std::cout << "histCol constructer " << towerId << std::endl;

  char name[] = "roccT00X17w3t4";
  char cvw[] = "XY";
  for( int layer=0; layer!=g_nLayer; layer++)
    for( int view=0; view!=g_nView; view++){
      for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
	for( int tDiv = 0; tDiv != g_nTime; tDiv++){
	  sprintf(name,"occT%d%c%dw%dt%d", 
		  towerId, cvw[view], layer, iWafer, tDiv);
	  nHits[layer][view][iWafer][tDiv] 
	    = new TH1F(name, name, g_nStrip, 0, g_nStrip);
	}
    }
}

void histCol::saveHists(){
  for( int layer=0; layer!=g_nLayer; layer++)
    for( int view=0; view!=g_nView; view++){
      for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
	for( int tDiv = 0; tDiv != g_nTime; tDiv++)
	  nHits[layer][view][iWafer][tDiv]->Write(0, TObject::kOverwrite);
    }
}
*/

//
// totCalib implementation 
//
totCalib::totCalib( const std::string jobXml, const std::string defJob ): 
  m_reconFile(0), m_reconTree(0), 
  m_reconEvent(0), m_digiFile(0), m_digiTree(0),
  m_digiEvent(0), m_rootFile(0)
{

  // get version number from CVS string
  std::string tag = "$Name:  $";
  int i = tag.find( " " );
  tag.assign( tag, i+1, tag.size() );
  i = tag.find( " " );
  tag.assign( tag, 0, i ) ;
  m_tag = tag;

  std::string version = "$Revision: 1.42 $";
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
  m_totAngleCF = 5.48E-1;
  m_peakMIP = 4.92;
  m_RSigma = 4.0;
  m_GFrac = 0.78;
  m_maxDirZ = -0.85;

  for(int tower = 0; tower != g_nTower; ++tower)
    m_towerPtr[tower] = -1;
  //
  // parse job options xml file
  //
  if( !readJobOptions( jobXml, defJob ) ) m_nEvents = -1;
  
  //
  std::cout << "totAngleCF: " << m_totAngleCF << ", peakMIP: " << m_peakMIP
	    << ", RSigma: " << m_RSigma
	    << ", GFrac: " << m_GFrac 
	    << ", maxDirZ: " << m_maxDirZ 
	    << std::endl;
  m_log << "totAngleCF: " << m_totAngleCF << ", peakMIP: " << m_peakMIP
	<< ", RSigma: " << m_RSigma
	<< ", GFrac: " << m_GFrac << ", maxDirZ: " << m_maxDirZ 
	<< std::endl;
  
}


void totCalib::getTimeStamp(){

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
}


bool totCalib::readJobOptions( const std::string jobXml, const std::string defJob ){
  
#ifdef OLD_RECON
  using namespace xml;
#else
  using namespace xmlBase;
#endif
  std::cout << "Open xml file: " << jobXml << std::endl;
  //
  // start parsing xml file.
  //
  XmlParser* parser = new XmlParser( true );
  DOMDocument* doc = 0;
  try {
    doc = parser->parse( jobXml.c_str() );
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    delete parser;
    return false;
  }
  if (doc != 0) {  // successful
    //std::cout << "Document successfully parsed" << std::endl;
    
    // top level element
    DOMElement* topElt = doc->getDocumentElement();
    //
    // get name of jobOption
    //
    std::string jobOption = defJob;
    if( jobOption == "None" ){
      try {
	jobOption = Dom::getAttribute(topElt, "default");
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
    }
    //
    // get default output directory
    //
    std::vector<DOMElement*> defOutList;
    Dom::getChildrenByTagName( topElt, "output", defOutList );

    //
    // get tot related parameters
    //
    std::vector<DOMElement*> totParamList;
    Dom::getChildrenByTagName( topElt, "totParam", totParamList );
    if( totParamList.size() > 0 ){
      DOMElement* totElt = totParamList.back();
      std::string value;
      try{
	value = Dom::getAttribute(totElt, "totAngleCF" );
	if( value.size() > 0 ) m_totAngleCF = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "peakMIP" );
	if( value.size() > 0 ) m_peakMIP = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "RSigma" );
	if( value.size() > 0 ) m_RSigma = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "GFrac" );
	if( value.size() > 0 ) m_GFrac = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "maxDirZ" );
	if( value.size() > 0 ) m_maxDirZ = atof( value.c_str() );
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
    }    

    //
    // search for jobOption
    //
    std::vector<DOMElement*> jobList;
    Dom::getChildrenByTagName( topElt, "jobOption", jobList );
    int len = jobList.size();
    if( len < 1 ){
      std::cout << "ERROR: no jobOption tag found in the jobOption xml file." << std::endl;
      return false;
    }
    
    bool jobFound = false;
    for(int i=0; i<len; i++){ //each jobOption loop
      DOMElement* jobElt = jobList[i];
      std::string name, type, mode="def";
      try{
	name = Dom::getAttribute(jobElt, "name");
	type = Dom::getAttribute(jobElt, "type");
	mode = Dom::getAttribute(jobElt, "mode");
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
      
      if( name != jobOption ) continue;
      
      //
      // jobOption found
      //
      jobFound = true;
      std::cout << "JobOption: " << jobOption << std::endl;

      // analysis type
      if( type == "badStrips" ){
	m_badStrips = true;
	std::cout << type << ", bad strip analysis" << std::endl;
      }
      else{
	m_badStrips = false;
	if( type == "MIP check" ){
	  m_correctedTot = true;
	  std::cout << type << ", tot check" << std::endl;
	}
	else{
	  m_correctedTot = false;
	  std::cout << type << ", tot analysis" << std::endl;
	}
      }

      // analysis mode
      if( mode == "hist" ){
	m_histMode = true;
	std::cout << "histogram only mode" << std::endl;
      }
      else{
	m_histMode = false;
	std::cout << "normal analysis mode" << std::endl;
      }

      // initialize root histograms
      initHists();
      
      // output
      std::string outDir;
      std::vector<DOMElement*> outList;
      Dom::getChildrenByTagName( jobElt, "output", outList );
      DOMNode* outElt;
      if( outList.size() > 0 )
	outElt = outList.back();
      else if( defOutList.size() > 0 )
	outElt = defOutList.back();
      else{
	std::cout << "no output directory specified." << std::endl;
	return false;
      }

      try{
	outDir = Dom::getAttribute(outElt, "dir");
	m_dtdDir = Dom::getAttribute(outElt, "dtdDir");
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
      if( m_dtdDir.size() < 5 ) // dtd directory is not specified. use default.
	m_dtdDir = "$(CALIBUTILROOT)/xml/";
      int status = facilities::Util::expandEnvVar(&m_dtdDir);
      if(status==-1){
	std::cout << m_dtdDir << " not found!" << std::endl;
	return false;
      }
      std::cout << "DTD directory: " << m_dtdDir << std::endl;
      
      getTimeStamp();
      if( !setOutputFiles( outDir.c_str() ) ) return false; 
      
      //
      // loop data tag
      // 
      TChain* digiChain = new TChain("Digi");
      TChain* reconChain = new TChain("Recon");
      std::vector<DOMElement*> dataList;
      Dom::getChildrenByTagName( jobElt, "data", dataList );
      int numData = (int) dataList.size();
      if( numData == 0 ){
	std::cout << "no data element found." << std::endl;
	m_log << "no data element found." << std::endl;
	return false;
      }
      std::string top, raw, digi, recon, runids;
      std::vector<std::string> runIds;
      for(int idata=0; idata<numData; idata++){ //each xml loop
	DOMNode* dataElt = dataList[ idata ];
	try{
	  top = Dom::getAttribute(dataElt, "top");
	  raw = Dom::getAttribute(dataElt, "raw");
	  digi = Dom::getAttribute(dataElt, "digi");
	  recon = Dom::getAttribute(dataElt, "recon");
	  runids = Dom::getAttribute(dataElt, "runIds");
	  mode = Dom::getAttribute(dataElt, "mode");
	}
	catch (DomException ex) {
	  std::cout << "DomException:  " << ex.getMsg() << std::endl;
	  return false;
	}
	runIds.clear();
	parseRunIds( runIds, runids );
	for(std::vector<std::string>::const_iterator run = runIds.begin();
	    run != runIds.end(); ++run) {

	  int runid = atoi( (*run).c_str() );
	  if( runid < m_first_run ) m_first_run = runid;
	  if( runid > m_last_run ) m_last_run = runid;
	}

	if( !readRcReports( raw.c_str(), runIds ) ) return false;
	if( m_badStrips )
	  if( !readHotStrips( raw.c_str(), runIds ) ) return false;
	if( mode == "dummy" ){
	  std::cout << "dummy mode: no digi/recon file used." << std::endl;
	  m_log << "dummy mode: no digi/recon file used." << std::endl;
	}
	else if( !addToChain( top.c_str(), digi.c_str(), recon.c_str(), 
			 runIds, digiChain, reconChain ) ) return false;
      }  
      m_nEvents =  setInputRootFiles( digiChain, reconChain );
      if( m_nEvents < 0 ) return false;

      //
      // loop hist tag
      //
      std::vector<DOMElement*> histList;
      Dom::getChildrenByTagName( jobElt, "hist", histList );
      int numHist = histList.size();
      for(int ihist=0; ihist<numHist; ihist++){ //each hist loop
	std::string dir, files;
	DOMNode* histElt = histList[ ihist ];
	try{
	  dir = Dom::getAttribute(histElt, "dir");
	  files = Dom::getAttribute(histElt, "files");
	}
	catch (DomException ex) {
	  std::cout << "DomException:  " << ex.getMsg() << std::endl;
	  return false;
	}
	std::vector<std::string> hfiles;
	splitWords( hfiles, files );
	if( !readInputHistFiles( dir, hfiles ) ) return false;
      }
      
      //
      // loop xml tag
      //
      std::vector<DOMElement*> xmlList;
      Dom::getChildrenByTagName( jobElt, "xml", xmlList );
      int numXml = xmlList.size();
      if( numXml == 0 ){
	std::cout << "no xml element found." << std::endl;
	m_log << "no xml element found." << std::endl;
	return false;
      }
      std::string atype, dir;
      for(int ixml=0; ixml<numXml; ixml++){ //each xml loop
	DOMNode* xmlElt = xmlList[ ixml ];
	try{
	  atype = Dom::getAttribute(xmlElt, "type");
	  dir = Dom::getAttribute(xmlElt, "dir");
	  runids = Dom::getAttribute(xmlElt, "runIds");
	}
	catch (DomException ex) {
	  std::cout << "DomException:  " << ex.getMsg() << std::endl;
	  return false;
	}
	if( atype != type ){
	  std::cout << " igonore inconsistent xml type: " << atype << "<>" << type << std::endl;
	  m_log << "ignore inconsistent xml type: " << atype << "<>" << type << std::endl;
	  continue;
	}
	runIds.clear();
	parseRunIds( runIds, runids );
	if( !readInputXmlFiles( dir, runIds ) ) return false;
      }      
    }
    if( !jobFound ){
      std::cout << "Invalid jobOption specified: " << jobOption << std::endl;
      return false;
    }
  }
  return true;
}


void totCalib::parseRunIds( std::vector<std::string>& runIds, 
			    const std::string& line ){

  splitWords( runIds, line );

  int runid;
  for( UInt_t i=0; i!=runIds.size(); i++){
    runid = atoi( runIds[i].c_str() );
    if( runid < 100000000 || runid > 999999999 ){
      std::cout << "Invalid run#: " << runid << std::endl;
      m_log << "Invalid run#: " << runid << std::endl;
      exit( EXIT_FAILURE );
    }
  }
}

void totCalib::splitWords(  std::vector<std::string>& words, 
		    const std::string& line ) {

  std::string::size_type pos = 0;
  std::string word;

  for( ; ; ) {

    string::size_type i = line.find(' ', pos);
    if(i != string::npos) word = line.substr(pos, i-pos);
    else word = line.substr(pos); // end of line
    words.push_back( word );
    
    if(i == string::npos) break;
    pos = i + 1;
  }
}

void totCalib::initHists(){
  m_nTrackDist = new TH1F("nTrack", "nTrack", 10, 0, 10);
  m_maxHitDist = new TH1F("maxHit", "maxHit", g_nUniPlane, 0, g_nUniPlane);
  m_numClsDist = new TH1F("numCls", "# of cluster per layer", 10, 0, 10 );
  m_dirzDist = new TH1F("dirZ", "dirZ", 100, -1, 1);
  m_armsDist = new TH1F("arms", "arms", 100, -5, 5);
  m_lrec = new TH1F("lrec", "lrec", g_nUniPlane, 0, g_nUniPlane);
  m_ldigi = new TH1F("ldigi", "ldigi", g_nUniPlane, 0, g_nUniPlane);
  m_lcls = new TH1F("lcls", "lcls", g_nUniPlane, 0, g_nUniPlane);
  if( m_badStrips ){
    m_locc = new TH1F("locc", "locc", g_nUniPlane, 0, g_nUniPlane);
    m_leff = new TH1F("leff", "leff", g_nUniPlane, 0, g_nUniPlane);
    m_ltrk = new TH1F("ltrk", "ltrk", g_nUniPlane, 0, g_nUniPlane);
    m_dist = new TH1F("dist", "distance", 50, 0, 200);
    m_brmsDist[0] = new TH1F("brms0", "brms 0-2", 100, -5, 5);
    m_brmsDist[1] = new TH1F("brms1", "brms 3-5", 100, -5, 5);
    m_brmsDist[2] = new TH1F("brms2", "brms 6-8", 100, -5, 5);
    m_brmsDist[3] = new TH1F("brms3", "brms 9-11", 100, -5, 5);
    m_brmsDist[4] = new TH1F("brms4", "brms 12-14", 100, -5, 5);
    m_brmsDist[5] = new TH1F("brms5", "brms 15-17", 100, -5, 5);
    m_occDist = new TH1F("occDist", "occDist", 200, 0, 200);
    m_poissonDist = new TH1F("poissonDist", "poissonDist", 40, -20, 0);
    m_aPos[0] = new TH1F("apos0", "apos0", 100, -250, 250);
    m_aPos[1] = new TH1F("apos1", "apos1", 100, -250, 250);
    m_aPos[2] = new TH1F("apos2", "apos2", 100, -250, 250);
    m_aPos[3] = new TH1F("apos3", "apos3", 100, -250, 250);
    m_aPos[4] = new TH1F("apos4", "apos4", 100, -250, 250);
    //m_aPos[0] = new TH1F("apos0", "apos0", 100, -50, 50);
    //m_aPos[1] = new TH1F("apos1", "apos1", 100, -50, 50);
    //m_aPos[2] = new TH1F("apos2", "apos2", 100, -50, 50);
    //m_aPos[3] = new TH1F("apos3", "apos3", 100, -50, 50);
  }
  else{
    m_fracBatTot = new TH1F("fracBadTot", "fraction of bad TOT", 50, 0, 0.2 );
    m_fracErrDist = new TH1F("fracErrDist", "Peak error", 100, 0, 0.1);
    m_chisqDist = new TH1F("chisqDist", "TOT fit chisq/ndf", 60, 0, 3);
    m_chargeScale = new TH1F("chargeScale", "Charge Scale", 50, 0.5, 1.5);
    m_langauWidth = new TH1F("langauWidth", "Langau Width", 50, 0.2, 0.7);
    m_langauGSigma = new TH1F("langauGSigma", "Langau GSigma", 50, 0.0, 2.0);
    m_dirProfile = new TProfile("dirProfile", "cons(theta) profile", 10, -1, -0.5);
    m_chist[4] = new TH1F("chargeAll", "TOT charge distribution", nTotHistBin, 0, maxTot);
    m_chist[0] = new TH1F("charge0", "TOT charge distribution (-1<cos<-0.95)", nTotHistBin, 0, maxTot);
    m_chist[1] = new TH1F("charge1", "TOT charge distribution (-0.95<cos<-0.9)", nTotHistBin, 0, maxTot);
    m_chist[2] = new TH1F("charge2", "TOT charge distribution (-0.9<cos<-0.85)", nTotHistBin, 0, maxTot);
    m_chist[3] = new TH1F("charge3", "TOT charge distribution (-0.85<cos<-0.8)", nTotHistBin, 0, maxTot);

  }

#ifdef FULLHIST
  for(int tower = 0; tower != g_nTower; ++tower) {
    for( int unp = 0; unp != g_nUniPlane; ++unp ){
      layerId lid( unp );
      int layer = lid.layer;
      int view = lid.view;
      char vw = 'X';
      if( view != 0 ) vw = 'Y';
      if( !m_badStrips ){
	char name[] = "var00000";
	sprintf(name,"var%d%d%d", tower, layer, view);
	m_totStrip[tower][layer][view] = new TGraphErrors(g_nDiv);
	m_totStrip[tower][layer][view]->SetName(name);
	
	char temp[] = "varCorr00000";
	sprintf(temp,"varCorr%d%d%d", tower, layer, view);
	m_chargeStrip[tower][layer][view] = new TGraphErrors(g_nDiv);
	m_chargeStrip[tower][layer][view]->SetName(temp);
	for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
	  char name1[] = "totT00X17fe0004";
	  sprintf(name1,"totT%d%c%dfe%d", tower, vw, layer, iDiv);
	  m_totHist[tower][layer][view][iDiv] = new TH1F(name1, name1, 100, 0, 200);
	  //char name2[] = "chargeT00X00fe0000";
	  //sprintf(name2,"chargeT%d%c%dfe%d", tower, vw, layer, iDiv);
	  //m_chargeHist[tower][unp][iDiv] = new TH1F(name2, name2, nTotHistBin, 0, maxTot);
	}
      }
    }
  }
#endif	    
}


totCalib::~totCalib() 
{
  if(m_rootFile == 0) return;

  std::cout << "save histgrams" << std::endl;

  m_rootFile->cd();

  m_nTrackDist->Write(0, TObject::kOverwrite);
  m_maxHitDist->Write(0, TObject::kOverwrite);
  m_numClsDist->Write(0, TObject::kOverwrite);
  m_dirzDist->Write(0, TObject::kOverwrite);
  m_armsDist->Write(0, TObject::kOverwrite);
  m_lrec->Write(0, TObject::kOverwrite);
  m_ldigi->Write(0, TObject::kOverwrite);
  m_lcls->Write(0, TObject::kOverwrite);
  if( m_badStrips ){
    for( int i=0; i<g_nLayer/3; i++) 
      m_brmsDist[i]->Write(0, TObject::kOverwrite);
    m_locc->Write(0, TObject::kOverwrite);
    m_leff->Write(0, TObject::kOverwrite);
    m_ltrk->Write(0, TObject::kOverwrite);
    m_dist->Write(0, TObject::kOverwrite);
    m_occDist->Write(0, TObject::kOverwrite);
    m_poissonDist->Write(0, TObject::kOverwrite);
    m_aPos[0]->Write(0, TObject::kOverwrite);
    m_aPos[1]->Write(0, TObject::kOverwrite);
    m_aPos[2]->Write(0, TObject::kOverwrite);
    m_aPos[3]->Write(0, TObject::kOverwrite);
    m_aPos[4]->Write(0, TObject::kOverwrite);
  }
  else{
    m_fracBatTot->Write(0, TObject::kOverwrite);
    m_fracErrDist->Write(0, TObject::kOverwrite);
    m_chisqDist->Write(0, TObject::kOverwrite);
    m_chargeScale->Write(0, TObject::kOverwrite);
    m_langauWidth->Write(0, TObject::kOverwrite);
    m_langauGSigma->Write(0, TObject::kOverwrite);
    m_dirProfile->Write(0, TObject::kOverwrite);
    for( int i=0; i!=5; i++)
      m_chist[i]->Write(0, TObject::kOverwrite);
    for( UInt_t i=0; i!=m_chargeHist.size(); i++ )
      m_chargeHist[i]->Write(0, TObject::kOverwrite);
  }

  for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw)
    m_towerVar[tw].saveHists();

  if( !m_badStrips ){
    for(int tower = 0; tower != g_nTower; ++tower) {
      for( int unp = 0; unp != g_nUniPlane; ++unp ){
#ifdef FULLHIST
	m_totStrip[tower][layer][view]->Write(0, TObject::kOverwrite);
	m_chargeStrip[tower][layer][view]->Write(0, TObject::kOverwrite);
#endif
	for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
#ifdef FULLHIST
	  m_totHist[tower][layer][view][iDiv]->Write(0, TObject::kOverwrite);
#endif
	  //m_chargeHist[tower][unp][iDiv]->Write(0, TObject::kOverwrite);
	}
      }
      
    }
  }

  m_rootFile->Close();
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
  m_rootFile = new TFile( filename.c_str(), "RECREATE" );
  if( m_rootFile ){
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
    exit( EXIT_FAILURE );
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
    exit( EXIT_FAILURE );
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
    exit( EXIT_FAILURE );
  }

  return nEvents;
}


int totCalib::setInputRootFiles( const char* rootDir, const char* digiPrefix,
				 const char* reconPrefix, 
				 const std::vector<std::string> &runIds )
{
  TChain* digiChain = new TChain("Digi");
  TChain* reconChain = new TChain("Recon");
  
  if( ! addToChain( rootDir, digiPrefix, reconPrefix, 
		    runIds, digiChain, reconChain ) ) return -1;
  
  return setInputRootFiles( digiChain, reconChain );
  
}

bool totCalib::addToChain( const char* rootDir, const char* digiPrefix,
			   const char* reconPrefix, 
			   const std::vector<std::string> &runIds,
			   TChain* digiChain, TChain* reconChain ){
  
  std::string digiFile, reconFile;
  char fname[] = "135000933/v4r060302p8/calib-v1r0/grRoot/recon-v3r1p2_135000933_recon_RECON_100.root";

  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {

    int runid = atoi( (*run).c_str() );

    digiFile = rootDir;
    sprintf(fname,"/%d/%s_%d_digi_DIGI.root",
	    runid, digiPrefix, runid);
    digiFile += fname;
    if( ! checkFile( digiFile ) ){
      std::cout << "digi file does not exist: " << digiFile << std::endl;
      m_log << "digi file does not exist: " << digiFile << std::endl;
      return false;
    }
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
	if( split == 0 ) {
	  std::cout << "recon file does not exist: " << reconFile << std::endl;
	  m_log << "recon file does not exist: " << reconFile << std::endl;
	  return false;
	}
	else break;
      }
      std::cout << "open recon file: " << reconFile << endl;
      m_log << "recon file: " << reconFile << endl;
      reconChain->Add( reconFile.c_str() );
      split++;
    }

  }
  return true;
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
    exit( EXIT_FAILURE );
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
    if( ! checkFile( reportFile ) ){
      std::cout << "Invalid rcReport file path: " << reportFile << std::endl;
      m_log << "Invalid rcReport file path: " << reportFile << std::endl;
      return false;
    }

    std::cout << "open rcReport file: " << reportFile << endl;
    m_log << "rcReport file: " << reportFile << endl;
    if( !parseRcReport( reportFile.c_str() ) ) return false;
  }
  return true;
}

bool totCalib::readHotStrips( const char* reportDir, 
				 const std::vector<std::string>& runIds )
{
  int runid=-100;
  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {
    runid++;
    if( atoi( (*run).c_str() ) == runid ) continue; // continuous from the previous.
    for( UInt_t tw=0; tw!=m_towerVar.size(); tw++){
      std::string xmlFile = reportDir;
      xmlFile += "/" + *run;
      xmlFile += "/" + m_towerVar[tw].hwserial;
      xmlFile += "_HotStrips.xml";
      //std::cout << xmlFile << std::endl;
      if( ! checkFile( xmlFile ) ){
	xmlFile = reportDir;
	xmlFile += "/" + *run;
	xmlFile += "/TkrHotStrips_" + m_towerVar[tw].hwserial;
	xmlFile += ".xml";
	//std::cout << xmlFile << std::endl;
	if( ! checkFile( xmlFile ) ) continue;
      }
      if( ! readBadStripsXmlFile( xmlFile.c_str() ) ) return false;
      runid = atoi( (*run).c_str() );
    }
    if( runid != atoi((*run).c_str()) ){ // no hot strips xml file found.
      std::cout << "no hot strip xml files found in " << reportDir 
		<< "/" << *run << std::endl;
      m_log << "no hot strip xml files found in " << reportDir 
	    << "/" << *run << std::endl;
      exit( EXIT_FAILURE );
    }
    
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

  if (docrcReport == 0){ //unsuccessful
    std::cout <<  "Parsing FAILURE: " << reportFile << std::endl;
    return false;
  }
  else{
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
      unsigned int pos = serials.find( "GTEM" );
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
	if( m_towerPtr[ tower ] < 0 ){
	  m_towerPtr[ tower ] = m_towerVar.size();
	  //towerVar dummy( tower );
	  m_towerVar.push_back( towerVar( tower, m_badStrips ) );
	  m_towerVar.back().hwserial = tower_serial;
	  std::cout << "Tower " << towerId << ": " << m_towerPtr[ tower ] 
		    << " " << tower_serial << std::endl;
	  m_log << "Tower " << towerId << ": " << m_towerPtr[ tower ] 
		    << " " << tower_serial << std::endl;
	}
	else if( m_towerVar[ m_towerPtr[tower] ].hwserial != tower_serial ){
	  std::cout << "Inconsistent tower serial IDs for tower " << tower
		    << ": " << m_towerVar[ m_towerPtr[tower] ].hwserial 
		    << " " << tower_serial
		    << std::endl;
	  m_log << "Inconsistent tower serial IDs for tower " << tower
		<< ": " << m_towerVar[ m_towerPtr[tower] ].hwserial 
		<< " " << tower_serial
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
    if( runid == m_first_run ){ 
      getDate( values[1].c_str(), m_startTime );
      std::cout << "start time: " << m_startTime 
		<< ", run id: " << runid << std::endl;
      m_log << "start time: " << m_startTime 
	    << ", run id: " << runid << std::endl;
    }
    if( runid == m_last_run ){
      getDate( values[2].c_str(), m_stopTime );
      std::cout << "stop time: " << m_stopTime 
		<< ", run id: " << runid << std::endl;
      m_log << "stop time: " << m_stopTime 
	    << ", run id: " << runid << std::endl;
    }
    
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

void totCalib::analyze( int nEvents )
{

  if( m_nEvents < 0 ) return; // something went wrong.
  else if( m_nEvents > 0 ){
    if( nEvents > 0 ) m_nEvents = nEvents;
    std::cout << "# of Events to analyze: " << m_nEvents << std::endl;
    m_log << "# of Events to analyze: " << m_nEvents << std::endl;

    analyzeEvents();
  }
  else{
    std::cout << "no events to analyze, skip root analysis." << std::endl;
    m_log << "no events to analyze, skip root analysis." << std::endl;
  }

  if( m_badStrips ){
    findBadStrips();
    if( ! m_histMode ) fillBadStrips();
  }
  else{
    fitTot();
    if( ! m_histMode ) fillXml();//takuya
  }
}

void totCalib::analyzeEvents() 
{
  int mEvent = int( m_nEvents * 0.01 );
  if( mEvent < 100 ) mEvent = 100;
  time_t startTime, currentTime;
  time( &startTime );

  for(int iEvent = 0; iEvent != m_nEvents; ++iEvent) {    
    if( iEvent >= mEvent ){
      time( &currentTime );
      int elapsedTime = currentTime - startTime;
      if( elapsedTime <= 0 ) elapsedTime = 1;
      int rEvents = m_nEvents - iEvent;
      if( elapsedTime > 2 )
	std::cout << "# of events: " << iEvent << " (" << iEvent*101/m_nEvents 
		  << "%) in " << elapsedTime << " s, "
		  << iEvent/elapsedTime << " events/s, "
		  << rEvents << " events, "
		  << int(1.0*rEvents*elapsedTime/iEvent)
		  << " s to go" << std::endl;
      if( mEvent > m_nEvents*0.095 ) mEvent += int( m_nEvents * 0.1 );
      else mEvent += int( m_nEvents * 0.01 );
    }
    if( m_reconEvent ) m_reconEvent->Clear();
    if( m_digiEvent ) m_digiEvent->Clear();

    m_reconTree->GetEntry(iEvent);
    m_digiTree->GetEntry(iEvent);

    assert(m_reconEvent != 0);
    assert(m_digiEvent != 0);

    if(! passCut()) continue;

    getReconClusters();
    getDigiClusters();
    selectGoodClusters();

    if( m_badStrips ) fillOccupancy( iEvent*g_nTime/m_nEvents );
    else fillTot();
  }
  std::cout << "Data scan finished." << std::endl;
  time( &currentTime );
  //protection against crash when testing with small number of events
  if(startTime==currentTime) return;
  std::cout << "total # of events: " << m_nEvents 
		<< " in " << (currentTime-startTime) << " s, "
		<< m_nEvents/(currentTime-startTime) << " events/s"
		<< std::endl;
  m_log << "total # of events: " << m_nEvents 
	<< " in " << (currentTime-startTime) << " s, "
	<< m_nEvents/(currentTime-startTime) << " events/s"
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

layerId totCalib::getLayerId( Cluster* cluster ){

  int tower = cluster->getTowerId();
  int unp = cluster->getUniPlane();

  layerId lid( unp );
  lid.setTower ( tower );
  return lid;
}


layerId totCalib::getLayerId( const TkrCluster* cluster )
{
#ifdef OLD_RECON
  int planeId = cluster->getPlane();
  TkrCluster::view viewId = cluster->getView();
  int tower = cluster->getTower();
  int layer = g_nLayer - planeId - 1;
  int view = (viewId == TkrCluster::X) ? 0 : 1;
#else
  commonRootData::TkrId id = cluster->getTkrId();
  int tower = TowerId( id.getTowerX(), id.getTowerY() ).id();
  int view = id.getView();
  int layer = cluster->getLayer();
#endif
  layerId lid( layer, view, tower);
  return lid;
}


void totCalib::getDigiClusters()
{
#ifdef PRINT_DEBUG
  std::cout << "getDigiClusters start" << std::endl;
#endif
  //
  // clear cluster information.
  for( UInt_t tw=0; tw!=m_towerVar.size(); tw++)
    for( int unp=0; unp!=g_nUniPlane; unp++)
      m_towerVar[tw].digiClusters[unp].clear();
  
  // The full collection of TkrDigis for this event
  const TObjArray* tkrDigiCol = m_digiEvent->getTkrDigiCol();
  if (!tkrDigiCol) return;
  
  std::vector<int> strips;
  // Loop over all TkrDigis
  TIter tkrIter(tkrDigiCol);
  TkrDigi *tkrDigi = 0;
  while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) {
    // Identify the tower and layer
    Int_t tower = tkrDigi->getTower().id();
    Int_t tw = m_towerPtr[ tower ];
    Int_t totl = tkrDigi->getToT(0);
    Int_t toth = tkrDigi->getToT(1);
    Int_t lastRC0Strip = tkrDigi->getLastController0Strip();

    Int_t layer = tkrDigi->getBilayer();
  
    // Returns the orientation of the strips
    GlastAxis::axis viewId = tkrDigi->getView();
    int view = (viewId == GlastAxis::X) ? 0 : 1;

    layerId lid( layer, view );
    int uniPlane = lid.uniPlane;
    
    strips.clear();
    UInt_t numHits = tkrDigi->getNumHits();
    // Loop through collection of hit strips for this TkrDigi
    UInt_t ihit;
    for (ihit = 0; ihit < numHits; ihit++) {
      // Retrieve the strip number
      Int_t iStrip = tkrDigi->getStrip(ihit);
      strips.push_back( iStrip );
      m_towerVar[tw].dHits[uniPlane][iStrip]++;      
    }
    std::sort( strips.begin(), strips.end() ); // sort strip# for clustering.
    for( UInt_t i=0; i!=strips.size(); i++){
      bool newCls = false;
      if( m_towerVar[tw].digiClusters[uniPlane].empty() ) newCls = true;
      else if( !m_towerVar[tw].digiClusters[uniPlane].back().addStrip( strips[i] ) )
	newCls = true;
      
      if( newCls ){
	m_ldigi->Fill( uniPlane );
	if( strips[i] <= lastRC0Strip )
	  m_towerVar[tw].digiClusters[uniPlane].push_back( Cluster( strips[i], totl ) );
	else
	  m_towerVar[tw].digiClusters[uniPlane].push_back( Cluster( strips[i], toth ) );
      }
    }
  }
#ifdef PRINT_DEBUG
  std::cout << "getDigiClusters end" << std::endl;
#endif
}

void totCalib::getReconClusters()
{
#ifdef PRINT_DEBUG
  std::cout << "getReconClusters start" << std::endl;
#endif
  
  // initialize recon cluster info
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++){
    int twr = m_towerPtr[ m_towerVar[tw].towerId ];
    if( twr != int(tw) ) {
      std::cout << "Invalid tower id: " << twr << " != " << tw << std::endl;
      m_log << "Invalid tower id: " << twr << " != " << tw << std::endl;
      exit( EXIT_FAILURE );
    }
    for( int unp=0; unp<g_nUniPlane; unp++)
      m_towerVar[tw].reconClusters[unp] = 0;
  }
  m_trackTowerList.clear();
  
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
  assert(tkrRecon != 0);
  
  int lastTower = -1;
  int numRecCls = 0;
  
#ifdef OLD_RECON
  std::map<int, TkrCluster*> clsMap;
  TObjArray* siClusterCol = tkrRecon->getClusterCol();
  int noOfTkrClusters = siClusterCol->GetLast()+1;
  for(int i = 0; i != noOfTkrClusters; ++i) {
    TkrCluster* cluster = dynamic_cast<TkrCluster*>(siClusterCol->At(i));
    clsMap[cluster->getId()] = cluster;
  }
  
  TkrKalFitTrack* tkrTrack = m_track;
  int nHitPlane = tkrTrack->getNumHits();
  for(int iPlane = 0; iPlane != nHitPlane; ++iPlane) {
    const TkrHitPlane* plane = tkrTrack->getHitPlane(iPlane);
    std::map<int, TkrCluster*>::const_iterator itr = clsMap.find(plane->getIdHit());
    assert(itr != clsMap.end());
    TkrCluster* cluster = itr->second;
#else
  TkrTrack* tkrTrack = m_track;
  TIter trk1HitsItr(tkrTrack);
  TkrTrackHit* pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {    
    const TkrCluster* cluster = (pTrk1Hit->getClusterPtr());
    if(!cluster) continue;
#endif
    numRecCls++;
    layerId lid = getLayerId( cluster );
    //std::cout << lid.tower << " " << lid.uniPlane << " " << lid.layer << " " << lid.view << std::endl;
    int tw = m_towerPtr[ lid.tower ];
    m_towerVar[tw].reconClusters[lid.uniPlane] = cluster;
    if( lid.tower != lastTower ){
      lastTower = lid.tower;
      m_trackTowerList.push_back( lastTower );
      if( lid.view == 0 )
	m_towerVar[tw].center[1] = (cluster->getPosition()).Y();
      else
	m_towerVar[tw].center[0] = (cluster->getPosition()).X();
    }
    m_lrec->Fill( lid.uniPlane );
  }
#ifdef PRINT_DEBUG
  std::cout << "getReconClusters end" << std::endl;
#endif
}

void totCalib::selectGoodClusters(){
  
#ifdef PRINT_DEBUG
  std::cout << "selectGoodClusters start" << std::endl;
#endif

  bool display = false;
  m_clusters.clear();
  //
  // register new raw clusters if it is close to the track position
  //
  Cluster* cluster;
  for( UInt_t tw=0; tw<m_trackTowerList.size(); tw++ ){
    int tower = m_trackTowerList[tw]; // order of towers for a track
    int twr = m_towerPtr[tower];
    for (int unp=g_nUniPlane-1; unp>=0; unp--){
      layerId lid( unp ); 
      int layer = lid.layer;
      int view = lid.view;
      float zpos = posZ[view][layer];
      //std::cout << layer << " " << view << " " 
      //	<< m_towerVar[twr].digiClusters[unp].size() << std::endl;

      int numCls = 0;
      for( UInt_t i=0; i!= m_towerVar[twr].digiClusters[unp].size(); i++){
	cluster = &( m_towerVar[twr].digiClusters[unp].at(i) );
	
	// calculate position
	float pos = ( cluster->getLastStrip() + cluster->getFirstStrip() )/2;
	int ladder = int( pos * 4 / g_nStrip );
	pos = pos - g_nStrip * 0.5 + 0.5;
	pos = pos * stripPitch + ladderGap * (ladder - 1.5 ) + m_towerVar[twr].center[view];
	//std::cout <<  cluster->getLastStrip() << " " 
	//<< cluster->getFirstStrip() << " " << pos << std::endl;
	if( view == 0 )
	  cluster->setXYZ( pos, m_towerVar[twr].center[1], zpos );
	else
	  cluster->setXYZ( m_towerVar[twr].center[0], pos, zpos );
	cluster->setId( tower, unp );
	if( m_towerVar[twr].reconClusters[unp] )
	  cluster->setCorrectedTot( 5.0 * m_towerVar[twr].reconClusters[unp]->getMips() );
	
	// check if this cluster is close to the track position
	layerId tlid;

	// find closest recon clusters
	float dzmin=10000, dzmin2=10000, dz;
	int numSkip=0, umin, umin2, tl, tunp;
	const TkrCluster* tcls;
	for( int dl=1; dl<g_nLayer; dl++){
	  for( int dir=-1; dir<2; dir+=2){
	    tl = layer + dl*dir;
	    if( tl < g_nLayer && tl >=0 ){
	      tlid.setLayer( tl, view );
	      tunp = tlid.uniPlane;
	      tcls = m_towerVar[twr].reconClusters[tunp];
	      if( tcls ){
		dz = fabs( zpos - tcls->getPosition().Z() );
		if( dz < dzmin ){
		  dzmin2 = dzmin;
		  umin2 = umin;
		  dzmin = dz;
		  umin = tunp;
		}
		else if( dz < dzmin2 ){
		  dzmin2 = dz;
		  umin2 = tunp;	
		}
		else numSkip++;
	      }
	    }
	    if( numSkip > 1 ) break;
	  }
	  if( numSkip > 1 ) break;
	}

	if( dzmin2>1000 || dzmin>1000 ) continue;

	TVector3 pos1 = m_towerVar[twr].reconClusters[umin]->getPosition();
	TVector3 pos2 = m_towerVar[twr].reconClusters[umin2]->getPosition();

	float delta;
	if( lid.view == 0 ){
	  pos = pos1.X() + ( pos2.X()-pos1.X() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
	  delta = cluster->getPosition().X() - pos;
	}
	else{
	  pos = pos1.Y() + ( pos2.Y()-pos1.Y() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
	  delta = cluster->getPosition().Y() - pos;
	}
	
	if( display ) std::cout << tower << " " << lid.layer << " " << lid.view << ", " << delta << " " << " " << pos << " " << zpos;
	m_armsDist->Fill( delta );
	if( fabs(delta) > 2.0 ){
	  if( display ) std::cout << " **************" << std::endl;
	  continue;
	}
	if( display ) std::cout << std::endl;
	
	//
	// good cluster, multiple cluster per layer allowed.
	//
	numCls++;
	m_clusters.push_back( cluster ); 
	for(int iStrip = cluster->getFirstStrip(); 
	    iStrip != int(cluster->getLastStrip()+1); ++iStrip){
	  m_towerVar[twr].rHits[unp][iStrip]++;
	  m_lcls->Fill( unp );
	}
      }
      m_numClsDist->Fill( numCls );
      //std::cout << tower << " " << uniPlane << std::endl;
    }
  }

#ifdef PRINT_DEBUG
  std::cout << "selectGoodClusters end" << std::endl;
#endif

}

bool totCalib::closeToTrack( const TkrCluster* cluster, TkrCluster* clusters[g_nTower][g_nUniPlane] )
{
  layerId lid = getLayerId( cluster ), tlid;
  int tower = lid.tower;
  int layer = lid.layer;
  int view = lid.view;
  float zpos = cluster->getPosition().Z();

  bool display = false;
  //if( layer == 4 && view == 0 ) display = true;
  if( display ) std::cout << "X4 ";

  // find closest hits
  float dzmin=10000, dzmin2=10000, dz;
  int numSkip=0, umin, umin2, tl, tunp;
  TkrCluster* tcls;
  for( int dl=1; dl<g_nLayer; dl++){
    for( int dir=-1; dir<2; dir+=2){
      tl = layer + dl*dir;
      if( tl < g_nLayer && tl >=0 ){
	tlid.setLayer( tl, view );
	tunp = tlid.uniPlane;
	tcls = clusters[tower][tunp];
	if( tcls ){
	  dz = fabs( zpos - tcls->getPosition().Z() );
	  if( dz < dzmin ){
	    dzmin2 = dzmin;
	    umin2 = umin;
	    dzmin = dz;
	    umin = tunp;
	  }
	  else if( dz < dzmin2 ){
	    dzmin2 = dz;
	    umin2 = tunp;	
	  }
	  else numSkip++;
	}
      }
      if( numSkip > 1 ) break;
    }
    if( numSkip > 1 ) break;
  }

  if( display ) std::cout << numSkip << " " << lid.uniPlane << " " 
			  << umin << " " << dzmin << " " 
			  << umin2 << " " << dzmin2 << ", ";
  
  if( dzmin2>1000 || dzmin>1000 ){
    if( display ) std::cout << std::endl;
    return false;
  }

  TVector3 pos1 = clusters[tower][umin]->getPosition();
  TVector3 pos2 = clusters[tower][umin2]->getPosition();

  float delta, pos;
  if( lid.view == 0 ){
    pos = pos1.X() + ( pos2.X()-pos1.X() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
    delta = cluster->getPosition().X() - pos;
  }
  else{
    pos = pos1.Y() + ( pos2.Y()-pos1.Y() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
    delta = cluster->getPosition().Y() - pos;
  }
  
  if( display ) std::cout << tower << " " << lid.layer << " " << lid.view << ", " << delta << " " << " " << pos << " " << zpos;
  m_armsDist->Fill( delta );
  if( fabs(delta) > 2.0 ){
    if( display ) std::cout << " **************" << std::endl;
    return false;
  }
  if( display ) std::cout << std::endl;

  return true;

}

int totCalib::findTot(int towerId, int layerId, int view , int stripId)
{
#ifdef PRINT_DEBUG
  std::cout << "findTot start" << std::endl;
#endif

  if(stripId <= m_lastRC0Strip[towerId][layerId][view] )
    return m_tot[towerId][layerId][view][0];
  else
    return m_tot[towerId][layerId][view][1];

#ifdef PRINT_DEBUG
  std::cout << "findTot end" << std::endl;
#endif
}

void totCalib::fillTot() 
{
  
#ifdef PRINT_DEBUG
  std::cout << "fillTot start" << std::endl;
#endif

  for( unsigned int cls=0; cls<m_clusters.size(); cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int unp = lid.uniPlane;
    int tw = m_towerPtr[ tower ];

    // require only a single strip
    if(cluster->getSize() != 1) continue;
    
    int iStrip = cluster->getFirstStrip();
    
    int tot = cluster->getRawToT();
    if( tot == 0 ) continue;
    
    float charge;
    if( m_correctedTot ) charge = cluster->getCorrectedTot();
    else charge  = calcCharge( lid, iStrip, tot);
    if( charge < 0.0 ) continue;
    charge /= ( 1 + m_totAngleCF * (1+m_dir.z()) ); // emprial correction factor
    m_chist[4]->Fill( charge );
    int idirz = int( 20 + m_dir.z()*20 );
    if( idirz < 4 ) m_chist[idirz]->Fill( charge );
    m_dirProfile->Fill( m_dir.z(), m_dir.z() );

    int iDiv = iStrip * g_nDiv / g_nStrip;
    int ibin = int( charge * nTotHistBin / maxTot );

    if( ibin < nTotHistBin && ibin >=0 )
      m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin]++;
    //m_chargeHist[tower][unp][iStrip/nStripPerGroup]->Fill(charge*(-m_dir.z()));
#ifdef FULLHIST
    m_totHist[tower][layer][view][iStrip/nStripPerGroup]->Fill(tot*(-m_dir.z())); 
#endif
  }
  
#ifdef PRINT_DEBUG
  std::cout << "fillTot end" << std::endl;
#endif
}

bool totCalib::passCut() 
{
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon(); 
  assert(tkrRecon != 0);
  
  TObjArray* tracks = tkrRecon->getTrackCol();
  int numTracks = tracks->GetEntries();
  m_nTrackDist->Fill( numTracks );
  
  // select only 1 or 2 track event
  if( numTracks > 2) return false;
  
  // find a track with maximum number of hits.
  int maxHits = 0, nHits;
  for( int tk=0; tk!=numTracks; tk++){
    
#ifdef OLD_RECON
    TkrKalFitTrack* track = dynamic_cast<TkrKalFitTrack*>(tracks->At(tk));
#else
    TkrTrack* track = dynamic_cast<TkrTrack*>(tracks->At(tk));
#endif
    if(track) {
      nHits = track->getNumFitHits();
      if( nHits > maxHits ){
	maxHits = nHits;
	m_track = track;
	m_pos = track->getInitialPosition();
	m_dir = track->getInitialDirection();
      }
    }
  }
  if( maxHits == 0 ) return false;
  m_maxHitDist->Fill( maxHits );
  m_dirzDist->Fill( m_dir.Z() );
  float maxDirZ = m_maxDirZ;
  if( m_badStrips ) maxDirZ = -0.7;
  if( m_dir.Z() > maxDirZ ) return false;
  
  return true;
}

void totCalib::fitTot()
{  
  // define Gaussian convolved Laudau function.
  TF1 *ffit = new TF1( "langau2", langau2fun, 0, 30, 6 );
  ffit->SetParNames( "Width", "MP", "Area", "GSigma", "RSigma", "GFrac" );
  std::cout << "Start fit." << std::endl;
  m_chargeHist.clear();
  
  const float meanChargeScale = 1.1, rangeChargeScale=0.3;
  char cvw[] = "XY";
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
    int tower = m_towerVar[ tw ].towerId;
    for(int unp = 0; unp != g_nUniPlane; ++unp) {
      layerId lid( unp );
      int layer = lid.layer;
      int view = lid.view;
      std::cout << "Tower " << tower << ": ";
      std::cout << cvw[view] << layer << std::endl;
      for(int iDiv = 0; iDiv != g_nDiv; ++iDiv){
	m_towerVar[tw].tcVar[unp].chargeScale[iDiv] = 1.0;
	
	float area, ave, rms;
	Double_t *par, *error;
	float peak, errPeak, width, errWidth;
	  
#ifdef FULLHIST
	// fit uncorrected tot for each strip
	area = m_totHist[tower][layer][view][iDiv]->Integral();
	ave = m_totHist[tower][layer][view][iDiv]->GetMean();
	rms = m_totHist[tower][layer][view][iDiv]->GetRMS();
	//std::cout << area << " " << ave << " " << rms << std::endl;
	if( area<100 || ave==0.0 || rms==0.0 ){ 
	  std::cout << "Layer: " << cvw[view] << layer
		    << ", FE: " << iDiv << ", Entries: " << area
		    << ", Mean: " << ave << ", RMS: " << rms 
		    << " skipped." << std::endl;
	  m_log << "Layer: " << cvw[view] << layer
		<< ", FE: " << iDiv << ", Entries: " << area
		<< ", Mean: " << ave << ", RMS: " << rms 
		<< " skipped." << std::endl;
	  continue;
	}
	
	ffit->SetParLimits( 0, 0.0, rms );
	ffit->SetParLimits( 1, 0.0, ave*2 );
	ffit->SetParLimits( 2, 0.0, area*0.4 );
	ffit->SetParLimits( 3, 0.0, rms );
	ffit->SetParLimits( 4, 1.0, 10.0 );
	ffit->SetParLimits( 5, 0.0, 1.0 );
	ffit->SetRange( ave-1.25*rms, ave+2*rms );
	ffit->SetParameters( rms*0.2, ave*0.75, area*0.1, rms*0.4 );
	ffit->FixParameter( 4, m_RSigma );
	ffit->FixParameter( 5, m_GFrac );
	//m_totHist[layer][view][iDiv]->Fit( "langau2", "RBQ" );
	
	//0:width(scale) 1:peak 2:total area 3:width(sigma)
	par = ffit->GetParameters();
	error = ffit->GetParErrors();
	
	float pos = float(iDiv);	  
	float errPos = 0.;
	
	peak = float( *(par+1) );
	errPeak = float( *(error+1) );
	
	width = float( *(par+3) );
	errWidth = float( *(error+3) );
	
	m_totStrip[tower][layer][view]->SetPoint(iDiv, pos, peak);
	m_totStrip[tower][layer][view]->SetPointError(iDiv, errPos, errPeak);
	
	/*m_log << "Uncorrected tot " << layer << ' ' << view << ' ' << pos
	  << ' ' << peak << ' ' << errPeak << ' ' << width << ' '
	  << errWidth << endl; */
#endif	  
	// fit charge for each strip
	char name[] = "chargeT00X00fe0000";
	sprintf(name,"chargeT%d%c%dfe%d", tower, cvw[view], layer, iDiv);
	TH1F* chargeHist = new TH1F(name, name, nTotHistBin, 0, maxTot);
	float binWidth = maxTot / nTotHistBin;
	for( int ibin=0; ibin!=nTotHistBin; ibin++)
	  chargeHist->Fill( (ibin+0.5)*binWidth, 
			    m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin] );
	m_chargeHist.push_back( chargeHist );

	area = chargeHist->Integral();
	ave = chargeHist->GetMean();
	rms = chargeHist->GetRMS();
	//std::cout << area << " " << ave << " " << rms << std::endl;
	if( area<200 || ave==0.0 || rms==0.0 ){ 
	  m_log << "T" << tower << " " << cvw[view] << layer
		<< " " << iDiv << ", Entries: " << area
		<< ", Mean: " << ave << ", RMS: " << rms 
		<< " skipped." << std::endl;
	  continue;
	}
	
	//
	// check the presence of bad TOT values
	// try not to include these TOT in the fit.
	//
	int bin = int(ave*0.5/binWidth) + 1;
	float fracBadTot = chargeHist->Integral(1,bin)
	  / chargeHist->Integral();
	m_fracBatTot->Fill( fracBadTot );

	float lowLim = ave - 1.25 * rms;
	if( fracBadTot > 0.05 && lowLim < ave*0.5 ){
	  lowLim = ave*0.5;
	  std::cout << "WARNING, large bad TOT fraction: " 
		    << fracBadTot << ", T" << tower << " " << cvw[view] 
		    << layer << " " << iDiv << std::endl;
	  m_log << "WARNING, large bad TOT fraction: " 
		<< fracBadTot << ", T" << tower << " " << cvw[view] << layer 
		<< " " << iDiv << std::endl;
	}

	ffit->SetParLimits( 0, 0.0, rms );
	ffit->SetParLimits( 1, 0.0, ave*2 );
	ffit->SetParLimits( 2, 0.0, area*0.4 );
	ffit->SetParLimits( 3, 0.0, rms );
	ffit->SetRange( lowLim, ave+2*rms );
	ffit->SetParameters( rms*0.2, ave*0.75, area*0.1, rms*0.4 );
	ffit->FixParameter( 4, m_RSigma );
	ffit->FixParameter( 5, m_GFrac );
	chargeHist->Fit( "langau2", "RBQ" );
	
	//0:width(scale) 1:peak 2:total area 3:width(sigma)
	par = ffit->GetParameters();
	error = ffit->GetParErrors();
	//par = (m_chargeHist[layer][view][iDiv]->GetFunction("landau"))->GetParameters();
	//error = (m_chargeHist[layer][view][iDiv]->GetFunction("landau"))->GetParErrors();
	
	peak = float( *(par+1) );
	errPeak = float( *(error+1) );
	if( peak > 0 ) m_fracErrDist->Fill( errPeak*sqrt(area/1000)/peak );
	
	width = float( *(par+3) );
	errWidth = float( *(error+3) );

	float chisq = ffit->GetChisquare();
	float ndf = ffit->GetNDF();
	if( ndf > 0 ) m_chisqDist->Fill( chisq/ndf );
#ifdef FULLHIST
	m_chargeStrip[tower][layer][view]->SetPoint(iDiv, pos, peak);
	m_chargeStrip[tower][layer][view]->SetPointError(iDiv, errPos, errPeak);
#endif	  
	if( peak > 0.0 ){
	  float chargeScale =  m_peakMIP / peak;
	  m_chargeScale->Fill( chargeScale );
	  m_langauWidth->Fill( *(par+0) ); //  width (scale)
	  m_langauGSigma->Fill( *(par+3) ); // width (sigma)
	  if( fabs(chargeScale-meanChargeScale) > rangeChargeScale ){
	    std::cout << "WARNING, Abnormal charge scale: " 
		      << chargeScale << ", T" << tower << " " << cvw[view] 
		      << layer << " " << iDiv << std::endl;
	    m_log << "WARNING, Abnormal charge scale: "
		  << chargeScale << ", T" << tower << " " << cvw[view] 
		  << layer << " " << iDiv << std::endl;
	    if( chargeScale > meanChargeScale ) 
	      chargeScale = meanChargeScale + rangeChargeScale;
	    if( chargeScale < meanChargeScale ) 
	      chargeScale = meanChargeScale - rangeChargeScale;
	  }
	  if( chisq/ndf > maxChisq ){ // large chisq/ndf
	    std::cout << "WARNING, large chisq/ndf: "
		      << chisq/ndf << ", T" << tower << " " << cvw[view] 
		      << layer << " " << iDiv << std::endl;
	    m_log << "WARNING, large chisq/ndf: "
		  << chisq/ndf << ", T" << tower << " " << cvw[view] << layer 
		  << " " << iDiv << std::endl;
	  }
	  // large peak fit error
	  if( errPeak*sqrt(area/1000)/peak > maxFracErr ){ 
	    std::cout << "WARNING, large peak fit error: "
		      << errPeak*sqrt(area/1000)/peak << ", T" << tower 
		      << " " << cvw[view] << layer << " " << iDiv << std::endl;
	    m_log << "WARNING, large peak fit error: "
		  << errPeak*sqrt(area/1000)/peak << ", T" << tower << " " 
		  << cvw[view] << layer << " " << iDiv << std::endl;
	  }
	  m_towerVar[tw].tcVar[unp].chargeScale[iDiv] = chargeScale;
	}
	
	m_log << "Fit T" << tower << " " << cvw[view] << layer << " " 
	      << iDiv << ' ';
	m_log.precision(3);
	m_log << area << ' ' << ave << ' ' << rms << ", " << *(par+0) << ' '
	      << *(par+1) << ' ' << *(par+2) << ' ' << *(par+3) << ", "
	      << errPeak/peak << " " << int(chisq+0.5) << "/" << ndf
	      << std::endl;
	
      }
    }
  }

  for( int i=0; i!=5; i++){
    float area = m_chist[i]->Integral();
    float ave = m_chist[i]->GetMean();
    float rms = m_chist[i]->GetRMS();
    ffit->SetParLimits( 0, 0.0, rms );
    ffit->SetParLimits( 1, 0.0, ave*2 );
    ffit->SetParLimits( 2, 0.0, area*0.4 );
    ffit->SetParLimits( 3, 0.0, rms );
    ffit->SetRange( ave-1.25*rms, ave+2*rms );
    ffit->ReleaseParameter( 4 );
    ffit->ReleaseParameter( 5 );
    ffit->SetParameters( rms*0.5, ave*0.75, area*0.1, rms*0.4, m_RSigma, m_GFrac );
    m_chist[i]->Fit( "langau2", "RBQ" );
  }
}


bool totCalib::readInputHistFiles(const std::string dir, 
				 const std::vector<std::string>& files ){

  for(UInt_t i=0; i!=files.size(); ++i) {
    string path;
    path = dir + files[i];
    m_log << "Open " << path << std::endl;
    std::cout << "Open " << path << std::endl;
    TFile* hfile = new TFile( path.c_str() );
    if( ! hfile ){
      std::cout << "File open failure: " << path << std::endl;
      m_log << "File open failure: " << path << std::endl;
      return false;
    }

    //
    // read histograms
    //
    if( !readHists( hfile, i, files.size() ) ) return false;

    // close hist file
    hfile->Close();
    delete hfile;
  }
  return true;
}

bool totCalib::readHists( TFile* hfile, UInt_t iRoot, UInt_t nRoot ){

  m_nTrackDist->Add( (TH1F*)hfile->Get( "nTrack" ) );
  m_maxHitDist->Add( (TH1F*)hfile->Get( "maxHit" ) );
  m_numClsDist->Add( (TH1F*)hfile->Get( "numCls" ) );
  m_dirzDist->Add( (TH1F*)hfile->Get( "dirZ" ) );
  m_armsDist->Add( (TH1F*)hfile->Get( "arms" ) );
  m_lrec->Add( (TH1F*)hfile->Get( "lrec" ) );
  m_ldigi->Add( (TH1F*)hfile->Get( "ldigi" ) );
  m_lcls->Add( (TH1F*)hfile->Get( "lcls" ) );
  if( m_badStrips ){
    char hname[]="brms0";
    for( int i=0; i<g_nLayer/3; i++){ 
      sprintf( hname, "brms%d", i );
      m_brmsDist[i]->Add( (TH1F*)hfile->Get( hname ) );
    }
    m_locc->Add( (TH1F*)hfile->Get( "locc" ) );
    m_leff->Add( (TH1F*)hfile->Get( "leff" ) );
    m_ltrk->Add( (TH1F*)hfile->Get( "ltrk" ) );
    m_dist->Add( (TH1F*)hfile->Get( "dist" ) );
    m_occDist->Add( (TH1F*)hfile->Get( "occDist" ) );
    m_poissonDist->Add( (TH1F*)hfile->Get( "poissonDist" ) );
    m_aPos[0]->Add( (TH1F*)hfile->Get( "apos0" ) );
    m_aPos[1]->Add( (TH1F*)hfile->Get( "apos1" ) );
    m_aPos[2]->Add( (TH1F*)hfile->Get( "apos2" ) );
    m_aPos[3]->Add( (TH1F*)hfile->Get( "apos3" ) );
    m_aPos[4]->Add( (TH1F*)hfile->Get( "apos4" ) );

    for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw)
      m_towerVar[tw].readHists( hfile, iRoot, nRoot );
  }
  else{
    m_fracBatTot->Add( (TH1F*)hfile->Get( "fracBadTot" ) );
    m_fracErrDist->Add( (TH1F*)hfile->Get( "fracErrDist" ) );
    m_chisqDist->Add( (TH1F*)hfile->Get( "chisqDist" ) );
    m_chargeScale->Add( (TH1F*)hfile->Get( "chargeScale" ) );
    m_langauWidth->Add( (TH1F*)hfile->Get( "langauWidth" ) );
    m_langauGSigma->Add( (TH1F*)hfile->Get( "langauSigma" ) );
    m_dirProfile->Add( (TH1F*)hfile->Get( "dirProfile" ) );
    for( int i=4; i!=-1; i--){
      char hname[]="chargeAll";
      if( i<4 ) sprintf( hname, "charge%d", i );
      m_chist[i]->Add( (TH1F*)hfile->Get( hname ) );
    }
    char cvw[] = "XY";
    for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
      int tower = m_towerVar[ tw ].towerId;
      for(int unp = 0; unp != g_nUniPlane; ++unp) {
	layerId lid( unp );
	int layer = lid.layer;
	int view = lid.view;
	for(int iDiv = 0; iDiv != g_nDiv; ++iDiv){
	  char name[] = "chargeT00X00fe0000";
	  sprintf(name,"chargeT%d%c%dfe%d", tower, cvw[view], layer, iDiv);
	  TH1F* hist = (TH1F*)hfile->Get( name );
	  for( int ibin=0; ibin!=nTotHistBin; ibin++)
	    m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin] 
	      += (int)hist->GetBinContent( ibin+1 );
	}
      }
    }
  }

  return true;
}

bool totCalib::readInputXmlFiles(const std::string dir, 
				 const std::vector<std::string>& runIds ){

  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {
    if( m_badStrips ){
      if( !readBadStripsXmlFile( dir.c_str(), (*run) ) )
	return false;
    }
    else if( !readTotConvXmlFile( dir.c_str(), (*run).c_str() ) )
      return false;
  }
  return true;
}


bool totCalib::readBadStripsXmlFile(const char* path, 
				    const std::string runid ){
  bool hotStrips = true;
  string filename;
  filename = path;
  if( runid != "None" ){
    hotStrips = false;
    char fname[] = "/398000364/TkrNoiseAndGain_398000364_dead.xml";
    sprintf(fname,"/%s/TkrNoiseAndGain_%s_dead.xml", runid.c_str(), runid.c_str() );
    filename += fname;
  }

  if( ! checkFile( filename ) ){
    std::cout << "Invalid bad strips xml file path: " << filename << std::endl;
    m_log << "Invalid bad strips xml file path: " << filename << std::endl;
    return false;
  }

  std::cout << "Open xml file: " << filename << std::endl;
  m_log << "Bad strips xml file: " << filename << std::endl;
  
  int length;
  char *buffer;
  std::string line;
  ifstream is( filename.c_str() );
  //is.open ( jobXml, ios::binary );
  while( is ){
    getline( is, line );
    if( line.substr(0,2) == "]>" ) break;
  }
  length = is.tellg();
  // get length of file:
  is.seekg (0, std::ios::end);
  length = int( is.tellg() ) - length;
  is.seekg (0, std::ios::beg);
  // allocate memory:
  buffer = new char [length];
  // read data as a block:
  is.seekg ( -length, std::ios::end );
  is.read ( buffer, length );
  is.close();
  //std::cout << line.substr(0,10) << " " << line.size() << std::endl;
  ofstream os;
  filename = "/tmp/temp" + m_timeStamp + ".xml";
  os.open( filename.c_str() );
  os.write( buffer, length );
  os.close();
  delete buffer;

#ifdef OLD_RECON
  //typedef xml xmlBase;
  using namespace xml;
#else
  using namespace xmlBase;
#endif

  XmlParser* parser = new XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse( filename.c_str() );
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    delete parser;
    return false;
  }
  os.open( filename.c_str() );
  os << "dummy";
  os.close();

  if (doc != 0) {  // successful
    //std::cout << "Document successfully parsed" << std::endl;

    // look up generic attributes
    DOMElement* docElt = doc->getDocumentElement();
    DOMElement* attElt = Dom::findFirstChildByName(docElt,"generic");
    
    std::string runId, hwserial;
    try {
      runId  = Dom::getAttribute(attElt, "runId");
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

    int towerId=-1, twr=-1;
    std::vector<std::string> stripList;
    for( unsigned int tw=0; tw<m_towerVar.size(); tw++)
      if( m_towerVar[tw].hwserial == hwserial ){
	if( !hotStrips ) m_towerVar[tw].runid = runId;
	towerId = m_towerVar[tw].towerId;
	twr = tw;
      }
    std::cout << "tower: " << towerId << ", serial: " << hwserial 
	      << ", runid: " << runId
	      << std::endl;
    if( towerId < 0 ) return false;

    XMLCh* xmlchElt = XMLString::transcode("uniplane");
    DOMNodeList* conList = doc->getElementsByTagName(xmlchElt);
    int len = conList->getLength();   

    for(int i=0; i<len; i++){ //loop each uniplane entry.
      DOMNode* childNode = conList->item(i);
      int tray = Dom::getIntAttribute(childNode, "tray");
      std::string which = Dom::getAttribute(childNode, "which");
      int howBad = Dom::getIntAttribute(childNode, "howBad");
      if( howBad > 4 ) continue;
      if( hotStrips && howBad>1 ) continue;
      if( !hotStrips && howBad<2 ) continue;
      //std::cout << "(tray,which)=(" << tray << ", " << which << ") howBad=" 
      //		<< howBad << std::endl;
     
      //get first child element
      DOMElement* elder = Dom::getFirstChildElement(childNode);
      DOMElement* younger;

      layerId lid( tray, which, towerId );
      int unp = lid.uniPlane;
      if( unp >= g_nUniPlane || unp < 0 ){
	std::cout << "Invalid uniPlane id: " << unp << std::endl;
	m_log << "Invalid uniPlane id: " << unp << std::endl;
	continue;
      }

      while( elder ){
	std::string strips = Dom::getAttribute( elder, "strips" );
	//std::cout << strips << std::endl;
	stripList.clear();
	splitWords( stripList, strips );
	for( UInt_t ist=0; ist!=stripList.size(); ist++){
	  if( stripList[ist].size() == 0 ) continue;
	  int strip = atoi( stripList[ist].c_str() );
	  if( strip < 0 || strip >= g_nStrip ){
	    std::cout << towerId << " " << tray << " " << which 
		      << ", invalid strip id: " << strip << std::endl;
	    continue;
	  }
	  // howBad=1; noisy, howbad=2; dead, howBad=4; disconnected.
	  m_towerVar[twr].bsVar[unp].knownBadStrips[howBad/2].push_back( strip );
	}
	younger = Dom::getSiblingElement( elder );
	elder = younger;
      }

    }
  }
  else return false;
  return true;
}

bool totCalib::readTotConvXmlFile(const char* dir, const char* runid)
{
  string filename;
  filename = dir;
  char fname[] = "/398000364/TkrTotGain_398000364.xml";
  sprintf(fname,"/%s/TkrTotGain_%s.xml", runid, runid);

  filename += fname;
  if( ! checkFile( filename ) ){
    std::cout << "Invalid TOT xml file path: " << filename << std::endl;
    m_log << "Invalid TOT xml file path: " << filename << std::endl;
    return false;
  }

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
    //std::cout << "Document successfully parsed" << std::endl;

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
    for( unsigned int tw=0; tw<m_towerVar.size(); tw++)
      if( m_towerVar[tw].hwserial == hwserial ){
	m_towerVar[tw].runid = tot_runid;
	towerId = m_towerVar[tw].towerId;
      }
    std::cout << "tower: " << towerId << ", serial: " << hwserial 
	      << ", runid: " << tot_runid
	      << std::endl;
    if( towerId < 0 ) return false;

    // look up tower attributes
    DOMElement* defElt = Dom::findFirstChildByName(attElt,"default");
    float intercept, slope, quad;
    if( defElt ){ // default attribute exists.
      try {
	intercept  = Dom::getDoubleAttribute(defElt, "intercept");
	slope  = Dom::getDoubleAttribute(defElt, "slope");
	quad  = Dom::getDoubleAttribute(defElt, "quad");
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
      }
    }

    XMLCh* xmlchElt = XMLString::transcode("uniplane");
    DOMNodeList* conList = doc->getElementsByTagName(xmlchElt);
    int len = conList->getLength();   
    if( len != g_nLayer*g_nView ){
      if( defElt ){
	int tw = m_towerPtr[ towerId ];
	for( int unp=0; unp!=g_nUniPlane; unp++)
	  for( int stripId =0; stripId!=g_nStrip; stripId++){
	    m_towerVar[tw].tcVar[unp].totOffset[stripId] = intercept;
	    m_towerVar[tw].tcVar[unp].totGain[stripId] = slope;
	    m_towerVar[tw].tcVar[unp].totQuadra[stripId] = quad;
	  }
	return true;
      }
      else{
	std::cout << "ERROR: # of layers in xml is invalid, " << len << std::endl;
	m_log << "ERROR: # of layers in xml is invalid, " << len << std::endl;
	return false;
      }
    }

    for(int i=0; i<len; i++){//each layers loop
      DOMNode* childNode = conList->item(i);
      int tray = Dom::getIntAttribute(childNode, "tray");
      std::string which = Dom::getAttribute(childNode, "which");
      //std::cout << "(tray,which)=(" << tray << ", " << which << ") ";
     
      //get first child element
      DOMElement* elder = Dom::getFirstChildElement(childNode);
      DOMElement* younger;

      layerId lid( tray, which, towerId );
      if( lid.uniPlane >= g_nUniPlane || lid.uniPlane < 0 ){
	std::cout << "Invalid uniPlane id: " << lid.uniPlane << std::endl;
	m_log << "Invalid uniPlane id: " << lid.uniPlane << std::endl;
	continue;
      }
	
      int numStrip = 0;
      while( getParam( elder, lid ) ){
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

bool totCalib::getParam(const DOMElement* totElement, layerId lid ){  
#ifdef OLD_RECON
  //typdef xml xmlBase;
  using namespace xml;
#else
  using namespace xmlBase;
#endif

  int tw = m_towerPtr[ lid.tower ];
  int unp = lid.uniPlane;
  int layer = lid.layer;
  int view = lid.view;

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
    //cout << "finished (layer,view)=(" << layer << ", "<< view << ")" << endl;
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
  
  m_towerVar[tw].tcVar[unp].totOffset[stripId] = values[0];
  m_towerVar[tw].tcVar[unp].totGain[stripId] = values[1];
  m_towerVar[tw].tcVar[unp].totQuadra[stripId] = values[2];
  return true;
}


float totCalib::calcCharge( layerId lid, int iStrip, int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  int tw = m_towerPtr[ lid.tower ];
  int unp = lid.uniPlane;
  // TOT to charge conversion
  float charge = m_towerVar[tw].tcVar[unp].totOffset[iStrip] 
    + time*m_towerVar[tw].tcVar[unp].totGain[iStrip]
    + time*time*m_towerVar[tw].tcVar[unp].totQuadra[iStrip];
  
  return charge;
}


void totCalib::fillXml()//takuya
{

  std::string dtdFile = m_dtdDir + "tkrCalibration.dtd";
  std::ifstream dtd( dtdFile.c_str() );
  if( dtd ){
    std::cout << "Open dtd file: " << dtdFile << std::endl;
    m_log << "dtd file: " << dtdFile << std::endl;
  }
  else{
    std::cout << dtdFile << " cannot be opened." << std::endl;
    return;
  }
  
  std::string filename = m_outputDir;
  char fname[] = "/TkrFMX_TkrChargeScale_050131-161012.xml";
  
  std::ofstream latxml, fmxml;
  std::string tot_runid = m_towerVar[0].runid;
  if( m_towerVar.size() > 1 ){
    sprintf( fname, "/LAT_TkrChargeScale_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    filename += fname;
    
    latxml.open( filename.c_str() );
    if( latxml ){
      std::cout << "Open LAT charge scale xml file: " << filename << std::endl;
      m_log << "LAT charge scale xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }
    for( unsigned int i=1; i<m_towerVar.size(); i++)
      tot_runid += ':' + m_towerVar[i].runid;
    openChargeScaleXml( latxml, dtd, tot_runid );
    latxml.precision(3);
  }

  for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
    int tower = m_towerVar[tw].towerId;
    filename = m_outputDir;
    sprintf( fname, "/%s_TkrChargeScale_%s-%s.xml", 
	     m_towerVar[tw].hwserial.c_str(), 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    filename += fname;

    fmxml.open( filename.c_str() );
    if( fmxml ){
      std::cout << "Open charge scale xml file: " << filename << std::endl;
      m_log << "Charge scale xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }
    tot_runid = m_towerVar[tw].runid;
    openChargeScaleXml( fmxml, dtd, tot_runid );
    fmxml.precision(3);
    fillTowerChargeScales( fmxml, tower );
    if( m_towerVar.size() > 1 ) fillTowerChargeScales( latxml, tower );

    fmxml << "</chargeScale>" << endl;
    fmxml.close();
  }
  if( m_towerVar.size() > 1 ){
    latxml << "</chargeScale>" << endl;
    latxml.close();
  }
  dtd.close();
}


void totCalib::fillTowerChargeScales( std::ofstream &xmlFile, const int tower ){
  TowerId twrId(tower); 
  int tower_row = twrId.iy();
  int tower_col = twrId.ix();
  int tw = m_towerPtr[ tower ];
  std::string hwserial = m_towerVar[tw].hwserial;
  xmlFile << "  <tower row=\"" << tower_row << "\" col=\"" << tower_col 
	  << "\" hwserial=\"" << hwserial << "\">" << endl;
  char cvw[] = "XY";
  
  for( int unp = 0; unp < g_nUniPlane; unp++){
    layerId lid( unp );
    int layer = lid.layer;
    int view = lid.view;
    int tray = lid.tray;
    std::string which = lid.which;

    xmlFile << std::endl
	    << "   <!-- **** layer " << cvw[view]  << layer 
	    << " **** -->" << std::endl
	    << "   <uniplane tray=\"" << tray << "\" which=\""
	    << which << "\">" << std::endl;
    for(int iDiv = 0; iDiv != g_nDiv; ++iDiv) {
      xmlFile << "    <gtfeScale id=\"" << iDiv << "\" chargeScale=\"" 
	      <<  m_towerVar[tw].tcVar[unp].chargeScale[iDiv] << "\"/>" << endl;
    }
    xmlFile << "   </uniplane>" << endl; 
  }
  xmlFile     << "  </tower>" << endl;  
}


void  totCalib::openChargeScaleXml( std::ofstream &xmlFile, std::ifstream &dtd, const std::string tot_runid ){

  xmlFile << "<?xml version=\"1.0\" ?>" << std::endl
	 << "<!DOCTYPE chargeScale [" << std::endl;

  std::string line;
  dtd.clear(); // reset status flag. This has to be done before moving pointer.
  dtd.seekg( 0, std::ios::beg ); // go back to the BOF.
  while( dtd ){
    getline(dtd, line);
    xmlFile << line << std::endl;
  }

  xmlFile << "]>" << std::endl;


  xmlFile << "<chargeScale>" << endl
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
  
}
  

void totCalib::fillOccupancy( int tDiv ) 
{
#ifdef PRINT_DEBUG
  std::cout << "fillOccupancy start" << std::endl;
#endif
  
  //initialize container
  int nHits[g_nTower][g_nLayer][g_nView][g_nWafer+1];
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++){
    int tower = m_towerVar[tw].towerId;
    for( int layer=0; layer<g_nLayer; layer++)
      for( int view=0; view<g_nView; view++)
	for( int i=0; i<g_nWafer+1; i++) nHits[tower][layer][view][i] = 0;
  }
  
  //
  // first loop to register hits and get tower offset
  //
  int hitLayers[g_nLayer];
  for( int layer=0; layer!=g_nLayer; layer++) hitLayers[layer]=0;
  
  for( unsigned int cls=0; cls<m_clusters.size(); cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int view = lid.view;
    int layer = lid.layer;
    
    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip){
      nHits[tower][layer][view][iStrip/384]++;
      nHits[tower][layer][view][g_nWafer]++;
    }
    
    hitLayers[layer]++;
    
  }
  
  //
  // main loop to fill occupancy and track position
  //
  float pos, apos, tpos[g_nView], lpos, posz;
  float dist, dxz, dyz, dz, delta, deltax, deltay;
  float dirX=m_dir.X()/m_dir.Z(), dirY=m_dir.Y()/m_dir.Z(), 
    preX=m_pos.X(), preY=m_pos.Y(), preXZ=m_pos.Z(), preYZ=m_pos.Z();
  int aview, preLayer=g_nLayer;
  int lastTower=-1, nTowers, towers[2];

  int numCls = m_clusters.size();
  for( int cls=0; cls<numCls; cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int view = lid.view;
    int layer = lid.layer;
    int unp = lid.uniPlane;
    
    // fill track positions in all layer between previous and current layer
    // carefull for moving across towers.
    int elyr = layer;
    if( cls == numCls-1 ) elyr = 0;
    for( int lyr=preLayer-1; lyr>= elyr; lyr--){ 
      // layers where hits are expected.
      if( hitLayers[lyr] != 0
	  || ( lyr!=0 && lyr!=g_nLayer-1 
	       && hitLayers[lyr+1]>0 && hitLayers[lyr-1]>0 ) ){
	posz =  posZ[view][lyr];
	dxz = posz - preXZ;
	tpos[0] = preX + dirX*dxz;
	dyz = posz - preYZ;
	tpos[1] = preY + dirY*dyz;
	
	// check if track moves across towers.
	if( tower != lastTower ){
	  if( lastTower < 0 ){ 
	    lastTower = tower;
	    towers[0] = tower;
	    nTowers = 1;
	  }
	  else{
	    towers[0] = lastTower;
	    towers[1] = tower;
	    nTowers = 2;
	    lastTower = -1;
	  }
	}
	for( int tw=0; tw<nTowers; tw++){
	  int twr = towers[tw];
	  int vtw = m_towerPtr[twr];
	  for( int vw=0; vw!=g_nView; vw++){
	    lpos = tpos[vw] - m_towerVar[vtw].center[vw];
	    for( int iw=0; iw!=g_nWafer; ++iw){
	      float stp = ( lpos-ladderGap*(iw-1.5) ) / stripPitch 
		+ g_nStrip/2;
	      if( stp>=iw*g_nStrip/4 && stp < (iw+1)*g_nStrip/4 ){
		layerId lid( lyr, vw );
		m_towerVar[vtw].bsVar[lid.uniPlane].tHits[int(stp)]++;
		m_ltrk->Fill( lid.uniPlane );
	      }
	    }
	  }
	} // tower loop
      } // valid layer
    } // layer loop
    preLayer = layer;

    TVector3 position = cluster->getPosition();
    posz = position.Z();
    if( fabs(posz-posZ[view][layer]) > 0.01 )
      std::cout << "Incosistent z position: " << layer << " " << view << ", " 
		<< posZ[view][layer] << " != " << posz << std::endl;
    dyz = posz - preYZ;
    dxz = posz - preXZ;
    deltax = preX + dirX*dxz - position.X();
    deltay = preY + dirY*dyz - position.Y();
    
    if( view == 0 ){
      dz = dxz;
      delta = deltax;
    }
    else{
      dz = dyz;
      delta = deltay;
    }
    
    float dx = dirX*dz;
    float dy = dirY*dz;
    dist =sqrt( dz*dz+dx*dx+dy*dy );
    m_dist->Fill( dist );
    if( dist < 30 ) dist = 30;
    delta *= (35.0/dist);
    m_brmsDist[layer/3]->Fill( delta );
    //if( layer==4 && view==0 ) m_brmsDist[layer/3]->Fill( delta );
    
    // select good clusters
    if( fabs(delta) > 3.0  ) continue;
    
    if( view == 0 ){
      aview = 1;
      pos = deltax;
      apos = deltay;
      if( dxz > 10.0 ) dirX = ( position.X() - preX ) / dxz;
      preX = position.X();
      preXZ = position.Z();
    }
    else{
      aview = 0;
      pos = deltay;
      apos = deltax;
      if( dyz > 10.0 ) dirY = ( position.Y() - preY ) / dyz;
      preY = position.Y();
      preYZ = position.Z();
    }
    
    //std::cout << layer << " " << view << ", " << pos << " " << apos
    //      << std::endl;
    
    int twr = m_towerPtr[tower];
    m_locc->Fill( lid.uniPlane );
    for(int iStrip = cluster->getFirstStrip(); 
	iStrip != int(cluster->getLastStrip()+1); ++iStrip){
      m_towerVar[twr].bsVar[unp].lHits[iStrip]++;
      if( nHits[tower][layer][aview][g_nWafer] > 0 ){
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( nHits[tower][layer][aview][iw] > 0 ){
	    m_towerVar[twr].bsVar[unp].nHits[iStrip][iw][tDiv]++;
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
      }
      else{
	for( int iw=0; iw<g_nWafer; iw++ )
	  if( fabs( apos-89.5*(iw-1.5) ) < 44.4 ){
	    m_towerVar[twr].bsVar[unp].nHits[iStrip][iw][tDiv]++;
	    m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
	  }
      }
    }
  }
  
#ifdef PRINT_DEBUG
  std::cout << "fillOccupancy end" << std::endl;
#endif
}


void totCalib::findBadStrips()
{

  float sumThreshold, occThreshold;

  m_leff->Add( m_locc );
  m_leff->Divide( m_ldigi );

  poissonFunc poisson;
  int occupancy;
  float maxProbD = -1000.0;

  for( unsigned int tw=0; tw!=m_towerVar.size(); ++tw){
    //int tower = m_towerVar[ tw ];
    int tower = m_towerVar[ tw ].towerId; 
    m_log << "Tower: " << tower << ", ID: " << m_towerVar[tw].hwserial << std::endl;
    for( int uniPlane = g_nUniPlane-1; uniPlane >=0; uniPlane--){
      layerId lid( uniPlane );
      lid.setTower ( tower );
      int layer = lid.layer;
      int view = lid.view;
      int unp = uniPlane;

      float eff = m_leff->GetBinContent( uniPlane+1 );
      char vw = 'X';
      if( view != 0 ) vw = 'Y';
      int numDeadStrips = 0;
      float meanRatio = 0.0;

      //
      // initialize known bad strip flag
      //
      bool known[g_nStrip];
      for( int i=0; i!=g_nStrip; i++)
	known[i] = false;
      //
      // flag known bad strip from online.
      //
      for( int iBad=0; iBad!=3; iBad++){
	if( m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size() == 0 ) continue;
	for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size(); i++){
	  int ist = m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i];
	  if( known[ist] ) continue;
	  known[ist] = true;
	}
      }
      //
      // determine rate fudge factor for this layer
      //
      int numGood = 0;
      float ratio;
      for( int strip=0; strip!=g_nStrip; strip++){
	if( m_towerVar[tw].bsVar[uniPlane].lHits[strip] == 0 ) continue;
	if( m_towerVar[tw].bsVar[uniPlane].tHits[strip] == 0 ) continue;
	ratio = m_towerVar[tw].bsVar[uniPlane].lHits[strip] 
	  / (m_towerVar[tw].bsVar[uniPlane].tHits[strip]*eff) ;
	if( fabs(ratio-1.6) > 0.4 ) continue;
	numGood++;
	meanRatio += ratio;
      }
      meanRatio /= numGood;
      eff *= meanRatio;
      m_log << "Layer: " << vw << layer << " " << eff << " " << meanRatio;
      if( eff < minEff ) m_log << "*****";
      m_log << std::endl;
	
      //
      // main bad strips finder
      //
      numGood = 0;
      meanRatio = 0.0;
      for( int strip=0; strip!=g_nStrip; strip++){

	//if( m_lHits[uniPlane][strip] == 0 ){ // dead strips
	if( m_towerVar[tw].bsVar[uniPlane].lHits[strip] == 0 ){ // dead strips
	  if( !known[strip] ) m_log << "new ";
	  m_log << strip << ": dead." << std::endl;
	  //categorize as disconnected
	  m_towerVar[tw].bsVar[uniPlane].badStrips[2].push_back(strip); 
	  m_towerVar[tw].bsVar[uniPlane].badStrips[g_nBad-1].push_back(strip); 
	  continue;
	}

	int fst = strip - g_nMerge;
	int nst = 1 + 2*g_nMerge;
	if( fst < 0 ) fst = 0;
	int est = fst + nst;
	if( est > g_nStrip ){
	  est = g_nStrip;
	  fst = est - nst;
	}
        float meanOcc = 0;
	for( int ist=fst; ist!=est; ist++)
	  //meanOcc += m_tHits[uniPlane][ist];
	  meanOcc += m_towerVar[tw].bsVar[uniPlane].tHits[ist];
        meanOcc *= (eff/nst);
        sumThreshold = meanOcc - 5*sqrt(meanOcc);
        if( sumThreshold < 2.0 ) sumThreshold = 2.0;
        float expOcc = meanOcc / (g_nTime*g_nWafer);
        //float expOcc = meanOcc / (g_nWafer);
        occThreshold = expOcc - 2*sqrt(expOcc);
        if( occThreshold < 1.1 ) occThreshold = 1.1;

	//
	// check if mean occupancy is correct.
	//
	float sum = 0.0, total=0.0, num=0;
	for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	  for( int iTime = 0; iTime != g_nTime; ++iTime ){
	    occupancy = m_towerVar[tw].bsVar[uniPlane].nHits[strip][iWafer][iTime];
	    total += occupancy;
	    if( occupancy < occThreshold ) continue;
	    num++;
	    sum += occupancy;
	  }
	}
	//m_log << strip << " " << sum << " " << num << " " << total 
	//      << "/" << meanOcc << "=" << total/meanOcc << std::endl;
	if( num > g_nWafer*g_nTime*0.4 || total > 0.5*meanOcc ){ 
	  float mean;
	  if( num > 0 ) mean = sum / num;
	  else mean = total / (g_nTime*g_nWafer);
	  if( mean<expOcc && mean>0.6*expOcc ){ // expected occupancy too large
	    expOcc = mean;
	    meanOcc = expOcc * (g_nTime*g_nWafer);
	    sumThreshold = meanOcc - 5*sqrt(meanOcc);
	    if( sumThreshold < 2.0 ) sumThreshold = 2.0;
	    occThreshold = expOcc - 2*sqrt(expOcc);
	    if( occThreshold < 1.1 ) occThreshold = 1.1;
	  }
	}

	sum = 0.0;
	bool occFlag = false, probFlag=false;
	int occ[g_nWafer][g_nTime], occW[g_nWafer], occT[g_nTime];
	float prob[g_nWafer][g_nTime], probW[g_nWafer], probT[g_nTime], 
	  probSum=0.0;
	for( int iTime = 0; iTime != g_nTime; ++iTime ){ 
	  probT[iTime] = 0.0;
	  occT[iTime] = 0;
	}
	for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	  int jWafer = iWafer;
	  if(layer%2 == 0) jWafer = g_nWafer-1-iWafer;
	  occW[ iWafer ] = 0;
	  probW[iWafer] = 0.0;
	  for( int iTime = 0; iTime != g_nTime; ++iTime ){
	    //occupancy = m_nHits[uniPlane][strip][jWafer][iTime];
	    occupancy = m_towerVar[tw].bsVar[uniPlane].nHits[strip][jWafer][iTime];
	    occ[iWafer][iTime] = occupancy;
	    occW[iWafer] += occupancy;
	    occT[iTime] += occupancy;
	    float probability = 0;
	    if( occupancy < occThreshold ){
	      probability = poisson.getLogIntProb( (double)expOcc, occupancy );
	    }
	    prob[iWafer][iTime] = probability;
	    if( occupancy ==0 || probability < maxProb ){
	      probW[iWafer] += probability;
	      probT[iTime] += probability;
	    }
	    else probT[iTime] = 0.0;
	  }
	  sum += occW[iWafer];
	  probSum += probW[iWafer];
	  m_occDist->Fill( occW[iWafer]+0.1 );
	  if( occW[iWafer] < occThreshold ) occFlag = true;
	  if( probW[iWafer] < maxProbW ) probFlag = true;
	}
	
	if( probSum > maxProbSum ) probFlag = false;

	if( !known[strip] && !occFlag && !probFlag && sum > sumThreshold ){
	  meanRatio += sum/meanOcc;
	  numGood++;
	  continue;
	}

	if( !known[strip] ) m_log << "!";
	m_log << strip << ": " << sum << "/" << (int)sumThreshold 
	      << " @" << (int)probSum << ", ";
	int nBad = 0;
	int mBad = 0;
	for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	  m_log << " " << occW[iWafer];
	}
	m_log << ", " << (int)occThreshold << " w";
	int numLowProbW = 0, numZeroOccW=0;
	float sumProbW = 0.0;
	for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
	  m_log << " " << (int)-probW[iWafer];
	  if( probW[iWafer] < maxProbW ){
	    numLowProbW++;
	    sumProbW += probW[iWafer];
	  }
	  else{
	    numLowProbW = 0;
	    sumProbW = 0.0;
	  }
	  if( occW[iWafer] == 0 ) numZeroOccW++;
	  else numZeroOccW = 0;
	}

	m_log << " t";
	probFlag = false;
	float sumProbT = 0.0, contSum=0.0;
	for(int iTime = 0; iTime != g_nTime; ++iTime){
	  m_log << " " << (int)-probT[iTime];
	  sumProbT += probT[iTime];
	  if( probT[iTime] < maxProbT ) contSum += probT[iTime];
	  else{ // not continuous
	    if( contSum < maxContSum ) probFlag = true; // evaluate current sum
	    contSum = 0.0;
	  }
	}

	// partially or intermittently disconnected
	//if( sumProbW < maxProbSum || sumProbT < maxProbSum ){
	if( sumProbT < maxProbSum || probFlag || contSum < maxContSum ){
	  if( numLowProbW == g_nWafer) nBad = 4;
	  else if( numLowProbW == numZeroOccW ) nBad = 3;
	  else nBad = 5;
	  m_log << " =" << nBad;
	}
	else if( sum < sumThreshold ){
	  m_log << "xxx";
	}
	if( nBad == 4 ) mBad = 4;

	sum = occW[0];
	m_log << ",";
	bool flagBreak = false;
	int occBreak = 0;
	for(int iWafer = 1; iWafer != g_nWafer; ++iWafer){
	  int value = occW[iWafer];
	  sum += value;
	  float mean = sum / (iWafer+1);
	  float probability = 0;
	  if( value < mean ) probability = poisson.getLogProb( mean, value );
	  m_log << " " << (int)-probability;
	  m_poissonDist->Fill( probability );
	  if( probability < probThreshold ){ 
	    flagBreak = true;
	    occBreak += occW[iWafer];
	  }
	  if( flagBreak ) numDeadStrips++;
	}
	if( flagBreak ){
	  if( occBreak == 0 ) mBad = 3; // partial disconnected
	  else mBad = 4; // intermittent partial disconnected
	}
	else if( nBad == 5 ) numDeadStrips += numLowProbW; 
	else if( nBad == 4 ) numDeadStrips += g_nWafer;
	
	if( !known[strip] && nBad <= 0 ){
	  m_log << std::endl;
	  meanRatio += sum/meanOcc;
	  numGood++;
	}
	else{
	  if( nBad > 0 ){
	    m_towerVar[tw].bsVar[uniPlane].badStrips[nBad].push_back(strip);
	    m_towerVar[tw].bsVar[uniPlane].badStrips[g_nBad-1].push_back(strip); 
	    m_log << " *" << mBad << " " << std::endl;
	    if( probSum > maxProbD ) maxProbD = probSum;
	  }
	  else m_log << " *miss*" << std::endl;

	  if( known[strip] || mBad != nBad ){
	    fixedDisp( (int)meanOcc/(g_nWafer*g_nTime), false );
	    for(int it=0; it<g_nTime; ++it) fixedDisp( occT[it] );
	    m_log << ",     ";
	    for(int it=0; it<g_nTime; ++it) fixedDisp( (int)-probT[it] );
	    m_log << std::endl;
	    for( int iw=0; iw!=g_nWafer; ++iw){
	      fixedDisp( occW[iw], false );
	      for(int it=0; it<g_nTime; ++it) fixedDisp( occ[iw][it] );
	      m_log << ",  ";
	      fixedDisp( (int)-probW[iw] );
	      for(int it=0; it<g_nTime; ++it) fixedDisp( (int)-prob[iw][it] );
	      m_log << std::endl;
	    }
	  }
	}
      }

      meanRatio /= numGood;

      //
      // combine with known bad channels
      //
      combineBadChannels( lid );

      m_log << vw << layer << ", # of bad channels:";
      std::cout << vw << layer << ", # of bad channel:";
      for( int iBad=0; iBad<g_nBad; iBad++){
	m_log << " " << m_towerVar[tw].bsVar[uniPlane].badStrips[iBad].size();
	std::cout << " " << m_towerVar[tw].bsVar[uniPlane].badStrips[iBad].size();
	//m_log << " " << m_badStrips[uniPlane][iBad].size();
	//std::cout << " " << m_badStrips[uniPlane][iBad].size();
      }
      m_log << "; " << m_towerVar[tw].bsVar[uniPlane].knownBadStrips[2].size() 
	    << ", " << numDeadStrips 
	    << ", " << meanRatio << std::endl;
      std::cout << "; " << m_towerVar[tw].bsVar[uniPlane].knownBadStrips[2].size() 
		<< ", " << numDeadStrips 
		<< ", " << meanRatio << std::endl;
    }
  }
}


void totCalib::combineBadChannels( layerId lid ){

  int tw = m_towerPtr[lid.tower];
  int unp = lid.uniPlane;
  //
  // initialize known bad strip flag
  //
  bool known[g_nStrip], dead[g_nStrip], hot[g_nStrip], found[g_nStrip], debug=false;
  for( int i=0; i!=g_nStrip; i++){
    known[i] = false;
    dead[i] = false;
    hot[i] = false;
    found[i] = false;
  }
  //
  // register bad strips found in this job (iBad=g_nBad-1)
  //
  int ist;
  char cvw[] = "XY";
  if( m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].size() != 0 ){
    for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].size(); i++){
      ist = m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i];
      if( found[ist] ){
	std::cout << "Redundant bad strip: T" << lid.tower << " " 
		  << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	m_log << "Redundant bad strip: T" << lid.tower << " " 
	      << cvw[lid.view] << lid.layer << " " << ist << std::endl;
      }
      found[ist] = true;
    }
  }
  //
  // add disconnected strip from online.
  //
  for( int iBad=2; iBad!=3; iBad++){
    if( m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size() == 0 ) continue;
    for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size(); i++){
      ist = m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i];
      if( known[ist] ) continue;
      if( !found[ist] ){
	//m_towerVar[tw].bsVar[unp].badStrips[4].push_back( ist );
	m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].push_back( ist );
	std::cout << "Missing known bad strip: T" << lid.tower << " " 
		  << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	m_log << "Missing known bad strip: T" << lid.tower << " " 
	      << cvw[lid.view] << lid.layer << " " << ist << std::endl;
      }
      known[ist] = true;
    }
  }
  //
  // register hot and dead strips (iBad=0,1) 
  // and tag them to be removed from disconnected strips.
  //
  for( int iBad=0; iBad!=2; iBad++){
    debug = false;
    if( m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size() == 0 ) continue;
    for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size(); i++){
      ist = m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i];
      if( iBad==0 && hot[ist] ) continue; // avoid duplicates
      if( iBad==1 && dead[ist] ) continue; // avoid duplicates
      if( !found[ist] ){
	std::cout << "Missing known dead strip: T" << lid.tower << " " 
		  << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	m_log << "Missing known dead strip: T" << lid.tower << " " 
	      << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	debug = true;
      }
      if( iBad == 0 ) hot[ist] = true;
      else            dead[ist] = true;
      m_towerVar[tw].bsVar[unp].badStrips[iBad].push_back( ist );
    }
    //
    // print debug information when inconsistency is found.
    //
    if( debug ){
      std::cout << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
		<< lid.layer << ", " << iBad << ": ";
      m_log << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
		<< lid.layer << ", " << iBad << ": ";
      for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].knownBadStrips[iBad].size(); i++){
 	std::cout << " " << m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i] 
		  << " " << known[m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i]];
 	m_log << " " << m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i] 
	      << " " << known[m_towerVar[tw].bsVar[unp].knownBadStrips[iBad][i]];
      }
      std::cout << ", ";
      m_log << ", ";
      for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].size(); i++){
 	std::cout << " " << m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i]
		  << " " << found[m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i]];
	m_log << " " << m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i]
	      << " " << found[m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i]];
      }
      std::cout << std::endl;
      m_log << std::endl;
    }

    sort( m_towerVar[tw].bsVar[unp].badStrips[iBad].begin(),
	  m_towerVar[tw].bsVar[unp].badStrips[iBad].end() );
  }

  //
  // remove hot and dead strips from disconnected strips.
  //
  for( int iBad=2; iBad!=g_nBad; iBad++){
    int numMatch = 0;
    debug = false;
    for( UInt_t i=0; i<m_towerVar[tw].bsVar[unp].badStrips[iBad].size(); i++){
      ist = m_towerVar[tw].bsVar[unp].badStrips[iBad][i];
      if( known[ist] ) continue; // leave disconnected strips found online
      if( !found[ist] ){
	std::cout << "Inconsistent dead strip: T" << lid.tower << " " 
		  << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	m_log << "Inconsistent bad strip: T" << lid.tower << " " 
	      << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	m_towerVar[tw].bsVar[unp].knownBadStrips[g_nBad-1].push_back( ist );
	debug = true;
      }
      if( (iBad!=g_nBad-1 && (dead[ist] || hot[ist])) ){ 
	m_towerVar[tw].bsVar[unp].badStrips[iBad][i] = g_nStrip;
	numMatch++;
      }
    }
    //
    // print debug information when inconsistency is found.
    //
    if( debug ){
      std::cout << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
		<< lid.layer << ", " << iBad << ": ";
      m_log << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
	    << lid.layer << ", " << iBad << ": ";
      for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].badStrips[iBad].size(); i++){
	std::cout << " " << m_towerVar[tw].bsVar[unp].badStrips[iBad][i];
	m_log << " " << m_towerVar[tw].bsVar[unp].badStrips[iBad][i];
      }
      std::cout << std::endl;
      m_log << std::endl;
    }

    if( numMatch > 0 ){ // remove matched strips
      sort( m_towerVar[tw].bsVar[unp].badStrips[iBad].begin(),
	    m_towerVar[tw].bsVar[unp].badStrips[iBad].end() );
      for( int i=0; i!=numMatch; i++){
	if( m_towerVar[tw].bsVar[unp].badStrips[iBad].back() != g_nStrip ){
	  std::cout << "Inconsistent matched strip: T" << lid.tower << " " 
		    << cvw[lid.view] << lid.layer << " " << ist << std::endl;
	  m_log << "Inconsistent matched strip: T" << lid.tower << " " 
		<< cvw[lid.view] << lid.layer << " " << ist << std::endl;
	  std::cout << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
		    << lid.layer << ", " << iBad << " " << numMatch << ": ";
	  m_log << "DEBUG: T" << lid.tower << " " << cvw[lid.view] 
		<< lid.layer << ", " << iBad << " " << numMatch << ": ";
	  for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].badStrips[iBad].size(); i++){
	    std::cout << " " << m_towerVar[tw].bsVar[unp].badStrips[iBad][i];
	    m_log << " " << m_towerVar[tw].bsVar[unp].badStrips[iBad][i];
	  }
	  std::cout << std::endl;
	  m_log << std::endl;
	  exit( EXIT_FAILURE );
	}
	m_towerVar[tw].bsVar[unp].badStrips[iBad].pop_back();
      }
    }
  }
}


void totCalib::fixedDisp( int ival, bool space ){
  //m_log.width(3);
  if( ival < 10 ) m_log << "  " << ival;
  else if( ival < 100 ) m_log << " " << ival;
  else if( space ) m_log << " **";
  else m_log << ival;
}


void totCalib::fillBadStrips()
{
  
  std::string dtdFile = m_dtdDir + "badStrips.dtd";
  std::ifstream dtd( dtdFile.c_str() );
  if( dtd ){
    std::cout << "Open dtd file: " << dtdFile << std::endl;
    m_log << "dtd file: " << dtdFile << std::endl;
  }
  else{
    std::cout << dtdFile << " cannot be opened." << std::endl;
    return;
  }
  
  std::string filename = m_outputDir;
  char fname[] = "/TkrFMX_DeadStrips_050131-161012.xml";

  std::ofstream latxml, fmxml;
  if( m_towerVar.size() > 1 ){
    sprintf( fname, "/LAT_DeadStrips_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    filename += fname;

    latxml.open( filename.c_str() );
    if( latxml ){
      std::cout << "Open LAT bad strips xml file: " << filename << std::endl;
      m_log << "LAT bad strips xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }
    openBadStripsXml( latxml, dtd );
  }

  for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
    filename = m_outputDir;
    sprintf( fname, "/%s_DeadStrips_%s-%s.xml", 
	     m_towerVar[tw].hwserial.c_str(), 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );  
    filename += fname;

    fmxml.open( filename.c_str() );
    if( fmxml ){
      std::cout << "Open bad strips xml file: " << filename << std::endl;
      m_log << "Bad strips xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }

    openBadStripsXml( fmxml, dtd );
    fillTowerBadStrips( fmxml, tw, g_nBad ); // add all bad sttrips field
    if( m_towerVar.size() > 1 ) fillTowerBadStrips( latxml, tw );

    fmxml << "</badStrips>" << std::endl;
    fmxml.close();
  }

  if( m_towerVar.size() > 1 ){
    latxml << "</badStrips>" << std::endl;
    latxml.close();
  }
  dtd.close();
}
  
  
void totCalib::openBadStripsXml( std::ofstream &xmlFile, std::ifstream &dtd ){
  
  xmlFile << "<?xml version=\"1.0\" ?>" << std::endl
	  << "<!DOCTYPE badStrips [" << std::endl;
  
  std::string line;
  dtd.clear(); // reset status flag. This has to be done before moving pointer.
  dtd.seekg( 0, std::ios::beg ); // go back to the BOF.
  while( dtd ){
    getline(dtd, line);
    xmlFile << line << std::endl;
  }
  
  xmlFile << "]>" << std::endl;
  
  
  xmlFile << "<badStrips badType=\"dead\">" << std::endl
	  << "<!-- includes partial dead strips; "
	  << " intermitent and/or wrire bond broken " 
	  << "between SSD wafers -->" << std::endl;
  
  xmlFile << "  <generic calType=\"stripOccupancy\" creatorName=\"badStrips\""
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
  
}
 

void totCalib::fillTowerBadStrips( std::ofstream &xmlFile, const int tw, 
				   const int nBad ){
  
  //dead, disconnected, partial disconnected, intermittently disconnected, 
  //intermittently partial disconnexcted
  int howBad[g_nBad] = {1,2,4,12,20,28,128}; 
  std::string cBad[g_nBad] = {"hot","dead","disconnected",
			      "partially disconnected",
			      "intermittently disconnected",
			      "intermittently partially connected",
			      "all bad (hot/online only included)"};
  char cvw[] = "XY";
  
  int tower = m_towerVar[tw].towerId;
  TowerId twrId( tower ); 
  int tower_row = twrId.iy();
  int tower_col = twrId.ix();
  std::string hwserial = m_towerVar[tw].hwserial;
  xmlFile << "  <tower row=\"" << tower_row << "\" col=\"" << tower_col 
	  << "\" hwserial=\"" << hwserial << "\""
	  << " nOnbdCalib=\"false\" nOnbdTrig=\"false\""
	  << " nOnbdData=\"false\"" << ">" << std::endl;
  
  for( int uniPlane = g_nUniPlane-1; uniPlane >=0; uniPlane--){
    layerId lid( uniPlane );
    int layer = lid.layer;
    int view = lid.view;
    int tray = lid.tray;
    std::string which = lid.which;
      
    xmlFile << std::endl
	    << "    <!-- layer " << cvw[view] << layer << " -->" << std::endl;
      
    for( int iBad=0; iBad!=nBad; iBad++ ){ 
      int itr = m_towerVar[tw].bsVar[uniPlane].badStrips[iBad].size();
      xmlFile << "    <!-- # of " << cBad[iBad] << " strips: " << itr 
	      << " -->" << std::endl;
      if( iBad == 0 ) xmlFile << "    <!-- " << std::endl; // hot strip is not included.
      xmlFile << "    <uniplane tray=\"" << tray << "\" which=\""
	      << which << "\" nOnbdCalib=\"false\" nOnbdTrig=\"false\""
	      << " nOnbdData=\"false\" howBad=\"" << howBad[iBad] << "\"";
      
      if(itr){
	xmlFile << ">" << std::endl << "      <stripList strips=\"";
	for(int i=0;i!=itr;i++)
	  xmlFile << " " << m_towerVar[tw].bsVar[uniPlane].badStrips[iBad][i];
	xmlFile << "\"/>" << std::endl 
		<< "    </uniplane>" << std::endl;
      }
      else xmlFile << "/>" << std::endl;
      if( iBad == 0 ) xmlFile << "    --> " << std::endl; // hot strip is not included.
      
    }
  }
  xmlFile << "  </tower>" << std::endl;
  
}

//
// poisson class implementation
//
poissonFunc::poissonFunc(){
  double fac=1.0; //prepare factorial upto 199
  float lfac=log(1.0); //prepare factorial upto 199
  factorial[0] = fac;
  logFactorial[0] = fac;

  for(int k=1;k!=200;k++){
    fac*=k;
    factorial[k] = fac;
    lfac += log(k);
    logFactorial[k] = lfac;
  }
}

double poissonFunc::getProb( double mean, int value ){
  double prob;
  if( mean > 100 || value > 100 ){
    float sigma = sqrt( mean );
    float delta = value - mean;
    prob = exp( delta*delta/(2*sigma) );
  }
  else{
    prob = pow( mean, (value-mean) ) 
      * ( factorial[(int)(mean+0.5)] / factorial[value] );
  }

  return prob;
}


float poissonFunc::getLogIntProb( double mean, int value ){
  int imean = (int)(mean+0.5);
  if( imean == 0 || value > imean ) return 0.0;
  // Poisson probability normalized by peak probability
  else if( value < 200 && imean < 200 ){
    double prob=0;
    for( int i=0; i!=value+1; i++ )
      prob += pow(mean,i)/factorial[i];
    return log10( prob ) - mean*log10( exp(1) );
  }
  // Gaussian probability
  else return - (value-mean)*(value-mean) / (2*mean);
}


float poissonFunc::getLogProb( float mean, int value ){
  int imean = (int)(mean+0.5);
  if( imean == 0 || value > imean ) return 0.0;
  // Poisson probability normalized by peak probability
  else if( value < 200 && imean < 200 ){
    return log(mean)*(value-imean) 
      + logFactorial[imean] - logFactorial[value];
    //float p_poisson = pow(mean,value)*exp(-mean)/factorial[value];
    //float norm      = pow(mean,imean)*exp(-mean)/factorial[imean];
    //return log( p_poisson/norm ); //normalize
  }
  // Gaussian probability
  else return - (value-mean)*(value-mean) / (2*mean);
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

//
// Two Gaussian convolved Lanbdau function
// it also employs variable integration step.
//
Double_t langau2fun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=Width2/Width of convoluted Gaussian function
   //par[5]=fraction of narrow Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  // variable interation step
  // integration step increases as points move to outside.
  Double_t np = 20.0;   // number of convolution steps in a region
  Double_t sc =  1.0;   // convolution extends to +-sc Gaussian sigmas
  Double_t sc2 =  4.0;  // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland, fgauss;
  Double_t sum = 0.0;
  Double_t xm, dx;
  Double_t step;
  Double_t sigma = par[3];
  Double_t sigma2 = par[4]*sigma;
  Double_t frac = par[5];
  Double_t width = par[0];
  if( sigma<=0.0 || width<=0.0 ) return 0.0;

  // MP shift correction
  mpc = par[1] - mpshift * width;
  if( mpc < 0.0 ) return 0.0;
  
  // Convolution integral parameters
  xm = x[0];
  Double_t dxmax = sc * sigma;
  step = dxmax * 2 / np;
  Double_t dxmax2 = sc2 * sigma2;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=0.5*step; dx<=dxmax2; dx+=step) {
    fgauss = step * ( frac * TMath::Gaus(0.0,dx,sigma) / sigma
                      + (1-frac) * TMath::Gaus(0.0,dx,sigma2) / sigma2 );
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;

    if( dx >= dxmax ){ // getting out of the curent region, change step
      dxmax *= 2;
      step *= 2;
    }
  }
  sum /= width;
  
  return ( par[2] * sum * invsq2pi );
}
