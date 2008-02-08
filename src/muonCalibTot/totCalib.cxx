
#include <cmath>
#include <ctime>
#include <cassert>

#include "totCalib.h"
#include "facilities/Util.h"
#include "commonRootData/idents/TowerId.h"

#include "calibTkrUtil/TkrHits.h"

using std::string;
using std::cout;
using std::endl;

XERCES_CPP_NAMESPACE_USE

//#define PRINT_DEBUG 1

float getTrunctedMean( std::vector<float> array, float fraction=0.9 ){

  if( array.size() == 0 ) return 0.0;

  sort( array.begin(), array.end() );
  int ib = array.size() * (1-fraction) * 0.5 + 0.5;
  int ie = array.size() - ib;

  if( ib>= ie ){
    ib = 0;
    ie = array.size();
  }

  float sum = 0;
  for( int index=ib; index<ie; index++) sum += array[index];
  return sum/(ie-ib);

}

ULong64_t axtoi(const char *hexStg) {
  int n = 0;         // position in string
  int m = 0;         // position in digit[] to shift
  int count = 16;         // loop index
  ULong64_t intValue = 0;  // integer value of hex string
  ULong64_t digit;      // hold values to convert
  while ( true ){
    if (hexStg[n]=='\0') break;
    n++;
  }
  if( n < 16 ) count = n;
  m = n - 1;
  n = 0;
  while(n < count) {
    if (hexStg[m] > 0x29 && hexStg[m] < 0x40 ) //if 0 to 9
      digit = hexStg[m] & 0x0f;            //convert to int
    else if (hexStg[m] >='a' && hexStg[m] <= 'f') //if a to f
      digit = (hexStg[m] & 0x0f) + 9;      //convert to int
    else if (hexStg[m] >='A' && hexStg[m] <= 'F') //if A to F
      digit = (hexStg[m] & 0x0f) + 9;      //convert to int
    else break;
    // digit is value of hex digit at position n
    // (n << 2) is the number of positions to shift
    // OR the bits into return value
    intValue = intValue | (digit << (n << 2));
    m--;   // adjust the position to set
    n++;   // next digit to process
  }
  return (intValue);
}


std::vector<int> bitFlagToList( ULong64_t bits ){

  std::vector<int> list;
  int bit = 0;
  while( bits ){
    ULong64_t digit = 1;
    digit = bits & (digit<<bit);
    if( digit ){
      bits -= digit;
      list.push_back( bit );
    }
    bit++;
  }
  return list;
}

//
// totCalib implementation 
//
totCalib::totCalib( const std::string jobXml, const std::string defJob ): 
  m_reconFile(0), m_reconTree(0), m_digiFile(0), m_digiTree(0)
{
  // get version number from CVS string
  std::string tag = "$Name:  $";
  int i = tag.find( " " );
  tag.assign( tag, i+1, tag.size() );
  i = tag.find( " " );
  tag.assign( tag, 0, i ) ;
  m_tag += ":" + tag;

  std::string version = "$Revision: 1.56 $";
  i = version.find( " " );
  version.assign( version, i+1, version.size() );
  i = version.find( " " );
  version.assign( version, 0, i ) ;
  m_version += ":" + version;
  std::cout << "Tag: " << m_tag << ", version: " << m_version << std::endl;

  m_first_run = 999999999;
  m_last_run = 0;
  m_startTime="01/20 2005, 22:52 GMT";
  m_stopTime="01/20 2005, 22:52 GMT";
  m_totThreshold = 1.177;
  m_totGain = 0.589;
  m_totQuad = 0.00490;

  //
  // parse job options xml file
  //
  if( !readJobOptions( jobXml, defJob ) ){
    m_nEvents = -1;
    std::cout << "Job Option file failure." << std::endl;
    exit( EXIT_FAILURE );
  }
  //
  std::cout << "totAngleCF: " << m_totAngleCF << ", peakMIP: " << m_peakMIP
	    << ", RSigma: " << m_RSigma
	    << ", GFrac: " << m_GFrac 
	    << ", maxDirZ: " << m_maxDirZ 
	    << ", maxTrackRMS: " << m_maxTrackRMS 
	    << ", maxDelta: " << m_maxDelta 
	    << ", totThreshold: " << m_totThreshold 
	    << ", totGain: " << m_totGain
	    << ", totQuad: " << m_totQuad 
	    << std::endl;
  m_log << "totAngleCF: " << m_totAngleCF << ", peakMIP: " << m_peakMIP
	<< ", RSigma: " << m_RSigma
	<< ", GFrac: " << m_GFrac << ", maxDirZ: " << m_maxDirZ 
	<< ", maxTrackRMS: " << m_maxTrackRMS 
	<< ", maxDelta: " << m_maxDelta 
	<< ", totThreshold: " << m_totThreshold 
	<< ", totGain: " << m_totGain
	<< ", totQuad: " << m_totQuad 
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
  
  using namespace xmlBase;
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
	value = Dom::getAttribute(totElt, "totThreshold" );
	if( value.size() > 0 ) m_totThreshold = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "totGain" );
	if( value.size() > 0 ) m_totGain = atof( value.c_str() );
	value = Dom::getAttribute(totElt, "totQuad" );
	if( value.size() > 0 ) m_totQuad = atof( value.c_str() );
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

      // TOT plots for all TFE
      if( !m_badStrips ) m_nDiv = g_nDiv;

      // initialize root histograms
      initCommonHists();
      if( m_badStrips ) initOccHists();
      initTotHists();

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
      std::string top, raw, digi, svac, recon, runids, dtype;
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
	  dtype = Dom::getAttribute(dataElt, "type");
	  m_nameType = Dom::getAttribute(dataElt, "nameType");
	  svac = Dom::getAttribute(dataElt, "svac");
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
	if( dtype == "mc" ){
	  char tower_serial[] = "TkrMC16";
	  for( int tower=0; tower!=g_nTower; tower++){
	    m_towerPtr[ tower ] = m_towerVar.size();
	    m_towerVar.push_back( towerVar( tower, m_badStrips ) );
	    sprintf( tower_serial, "TkrMC%d", tower );
	    m_towerVar.back().hwserial = tower_serial;
	    std::cout << "Tower " << tower << ": " << m_towerPtr[ tower ] 
		      << " " << tower_serial << std::endl;
	    m_log << "Tower " << tower << ": " << m_towerPtr[ tower ] 
		  << " " << tower_serial << std::endl;
	    if( !m_badStrips )
	      for( int unp=0; unp!=g_nUniPlane; unp++)
		for( int stripId =0; stripId!=g_nStrip; stripId++){
		  m_towerVar.back().tcVar[unp].totThreshold[stripId] = m_totThreshold;
		  m_towerVar.back().tcVar[unp].totGain[stripId] = m_totGain;
		  m_towerVar.back().tcVar[unp].totQuad[stripId] = m_totQuad;
		}
	  }
	}
	else{
	  if( !readRcReports( raw.c_str(), runIds ) ) return false;
	  if( m_badStrips )
	    if( !readHotStrips( raw.c_str(), runIds ) ) return false;
	}
	if( mode == "dummy" ){
	  std::cout << "dummy mode: no digi/recon file used." << std::endl;
	  m_log << "dummy mode: no digi/recon file used." << std::endl;
	}
	else if( mode == "svac" ){ // read TKR histograms from svac root files.
	  if( !readInputHistFiles( top, svac, runIds ) ) return false;
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
      if( dtype != "mc" && numXml == 0 ){
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
      //
      // loop text tag
      //
      std::string names;
      std::vector<DOMElement*> txtList;
      Dom::getChildrenByTagName( jobElt, "text", txtList );
      int numTxt = txtList.size();
      for(int itxt=0; itxt<numTxt; itxt++){ //each text loop
	DOMNode* txtElt = txtList[ itxt ];
	try{
	  atype = Dom::getAttribute(txtElt, "type");
	  dir = Dom::getAttribute(txtElt, "dir");
	  names = Dom::getAttribute(txtElt, "names");
	}
	catch (DomException ex) {
	  std::cout << "DomException:  " << ex.getMsg() << std::endl;
	  return false;
	}
	if( atype != "knownBadStrips" ){
	  std::cout << " igonore inconsistent txt type: " << atype << "<>" << type << std::endl;
	  m_log << "ignore inconsistent txt type: " << atype << "<>" << type << std::endl;
	  continue;
	}
	if( !readBadStripsTxtFiles( dir, names ) ) return false;
      }
      break; // stop after reading the first matching job option file
    }
    if( !jobFound ){
      std::cout << "Invalid jobOption specified: " << jobOption << std::endl;
      return false;
    }
  }
  delete parser;
  delete doc;
  return true;
}


void totCalib::parseRunIds( std::vector<std::string>& runIds, 
			    const std::string& line ){

  splitWords( runIds, line );

  int runid;
  for( UInt_t i=0; i!=runIds.size(); i++){
    runid = atoi( runIds[i].c_str() );
    if( runid < 0 || runid > 999999999 ){
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


totCalib::~totCalib() 
{
  if(m_rootFile == 0) return;

  m_rootFile->cd();
  saveAllHist( true, false );

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
  if( m_log.is_open() )
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

    digiFile = rootDir;
    if( m_nameType == "straight" ) sprintf(fname,"/%s_%s_digi.root",
					   digiPrefix, (*run).c_str());
    else sprintf(fname,"/%s/%s_%s_digi_DIGI.root",
		 (*run).c_str(), digiPrefix, (*run).c_str());
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
	if( m_nameType == "straight" ) sprintf(fname,"/%s_%s_recon.root",
					       reconPrefix, (*run).c_str());
	else sprintf(fname,"/%s/%s_%s_recon_RECON.root",
		     (*run).c_str(), reconPrefix, (*run).c_str());
      else
	if( m_nameType == "straight" ) 
	  sprintf(fname,"/%s_%s_recon_%d.root",
		  reconPrefix, (*run).c_str(), split);
	else sprintf(fname,"/%s/%s_%s_recon_RECON_%d.root",
		     (*run).c_str(), reconPrefix, (*run).c_str(), split);
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
    if( ! checkFile( reportFile ) ){ // check LICOS directopry structure
      reportFile = reportDir;
      reportFile += "/" + *run;
      reportFile += "/LICOS/rcReport.out";
    }
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
    if( abs(atoi( (*run).c_str() )-runid)<3 ){
      runid = atoi( (*run).c_str() );
      continue; // continuous from the previous.
    }
    bool found = true;
    if( atoi( (*run).c_str() ) < 100000000 ){ // LICOS 
      std::string xmlFile = reportDir;
      xmlFile += "/" + *run;
      xmlFile += "/LICOS/latc_intDefaultLatc_TFE.xml";
      std::cout << xmlFile << std::endl;
      found = checkFile( xmlFile );
      if( found )
	if( ! readLatcTfeXmlFile( xmlFile ) ) return false;
    }
    else{
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
	if( ! readBadStripsXmlFile( xmlFile, true ) ) return false;
      }
      found = true;
    }
    runid = atoi( (*run).c_str() );
    if( !found ) { // no hot strips xml file found.
      std::cout << "no hot strip xml files found in " << reportDir 
		<< "/" << *run << std::endl;
      m_log << "no hot strip xml files found in " << reportDir 
	    << "/" << *run << std::endl;
      //exit( EXIT_FAILURE );
    }
    
  }
  return true;
}


bool totCalib::parseRcReport( const char* reportFile )
{
  using namespace xmlBase;

  XmlParser* parsercReport = new XmlParser(true);
  DOMDocument* docrcReport = 0;
  try{
    docrcReport = parsercReport -> parse(reportFile);
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
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
    //keywords.push_back("SerialNos");

    for( unsigned int i=0; i<keywords.size(); i++){
      try {
	DOMElement* childElt 
	  = Dom::findFirstChildByName( rcElt, keywords[i].c_str() );
	values.push_back( Dom::getTextContent(childElt) );
	//std::cout << keywords[i] << ": " << values[i] << std::endl;
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
	return false;
      }
    }

    int runid = atoi( values[0].c_str() );
    std::string serials, towerId, tower_serial;
    if( runid < 100000000 ){ // LICOS, use known serial IDs.
      serials = "{'GTEM(0,)': {'tkr': 'TkrFMA'}, 'GTEM(1,)': {'tkr': 'TkrFM2'}, 'GTEM(2,)': {'tkr': 'TkrFM14'}, 'GTEM(3,)': {'tkr': 'TkrFM15'}, 'GTEM(4,)': {'tkr': 'TkrFMB'}, 'GTEM(5,)': {'tkr': 'TkrFM1'}, 'GTEM(6,)': {'tkr': 'TkrFM12'}, 'GTEM(7,)': {'tkr': 'TkrFM13'}, 'GTEM(8,)': {'tkr': 'TkrFM5'}, 'GTEM(9,)': {'tkr': 'TkrFM3'}, 'GTEM(10,)': {'tkr': 'TkrFM7'}, 'GTEM(11,)': {'tkr': 'TkrFM9'}, 'GTEM(12,)': {'tkr': 'TkrFM6'}, 'GTEM(13,)': {'tkr': 'TkrFM4'}, 'GTEM(14,)': {'tkr': 'TkrFM10'}, 'GTEM(15,)': {'tkr': 'TkrFM11'}}";
    }
    else{ // normal runs, get serials from rcReport
      serials = Dom::getTextContent( Dom::findFirstChildByName( rcElt, "SerialNos" ) );
    }

    while( true ){
      unsigned int pos = serials.find( "GTEM" );
      if( pos == string::npos ) break;
      serials.assign( serials, pos+5, serials.size() );
      pos = serials.find( "," );
      towerId.assign( serials, 0, pos );
      pos = serials.find( "tkr" );
      unsigned int pos2 = serials.find( "GTEM" );
      if( pos2 != string::npos && pos2 < pos ) continue;
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

    if( runid != m_first_run && runid != m_last_run ){
      delete docrcReport;
      delete parsercReport;
      return true;
    }
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
    
    delete docrcReport;
    delete parsercReport;
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
    if( nEvents > 0 && nEvents < m_nEvents ) m_nEvents = nEvents;
    std::cout << "# of Events to analyze: " << m_nEvents << std::endl;
    m_log << "# of Events to analyze: " << m_nEvents << std::endl;

    analyzeEvents();
  }
  else{
    std::cout << "no events to analyze, skip root analysis." << std::endl;
    m_log << "no events to analyze, skip root analysis." << std::endl;
  }

  fitTot();
  if( m_badStrips ){
    if( ! m_histMode ){ 
      findBadStrips();
      fillBadStrips();
    }
    //calculateEfficiency();
  }
  else if( !m_histMode ) fillXml();//takuya

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

    monitorTKR();

    if(! passCut()) continue;

    getReconClusters();
    getDigiClusters();
    selectGoodClusters();

    if( m_badStrips ) fillOccupancy( iEvent*g_nTime/m_nEvents );
    fillTot();
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



bool totCalib::readInputHistFiles(const std::string dir, 
				 const std::vector<std::string>& files ){

  for(UInt_t i=0; i!=files.size(); ++i) {
    string path;
    path = dir + '/' + files[i];
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

bool totCalib::readInputHistFiles( const std::string dir, 
				   const std::string prefix, 
				   const std::vector<std::string> &runIds ){
  
  std::string path;
  char fname[] = "/135005345/v05r0703p06/calib-v1r0/svacTuple/emRootv0r0/svacTuple-v03r04p04_135005345_svac_svac.root";

  for(UInt_t i=0; i!=runIds.size(); ++i) {
    path = dir;
    sprintf(fname,"/%s/%s_%s_svac_svac.root",
	    runIds[i].c_str(), prefix.c_str(), runIds[i].c_str() );
    path += fname;
    m_log << "Open " << path << std::endl;
    std::cout << "Open " << path << std::endl;
    TFile* hfile = new TFile( path.c_str() );
    if( ! hfile ){
      std::cout << "File open failure: " << path << std::endl;
      m_log << "File open failure: " << path << std::endl;
      return false;
    }
    hfile->cd( "TkrCalib" );
    //
    // read histograms
    //
    if( !readHists( hfile, i, runIds.size() ) ) return false;

    // close hist file
    hfile->Close();
    delete hfile;
  }
  return true;
}

template <class HIST> void addHIST( HIST* hist, TFile* hfile, char* name  ){
  //std::cout << name << hist << std::endl;
  HIST* th1 =  (HIST*)hfile->FindObjectAny( name );
  if( th1 && hist ){
    hist->Add( th1 );
    //std::cout << th1->GetEntries() << std::endl;
    //std::cout << hist->GetEntries() << std::endl;
  }
  else{
    std::cout << "Invalid Histogram: " << name << std::endl;
    //m_log << "Invalid Histogram: " << name << std::endl;
  }
  
};

bool totCalib::readHists( TFile* hfile, UInt_t iRoot, UInt_t nRoot ){

  addHIST( m_nTrackDist, hfile, "nTrack" );
  addHIST( m_maxHitDist, hfile, "maxHit" );
  addHIST( m_trkRMS, hfile, "trkRMS" );
  addHIST( m_rawTOT, hfile, "rawTOT" );
  addHIST( m_numHitGTRC, hfile, "numHitGTRC" );
  addHIST( m_trkRMS1TWR, hfile, "trkRMS1TWR" );
  addHIST( m_trkRMS2TWR, hfile, "trkRMS2TWR" ); 
  addHIST( m_rmsProf1TWR, hfile, "rmsProf1TWR" );
  addHIST( m_rmsProf2TWR, hfile, "rmsProf2TWR" );
  addHIST( m_tresProfX, hfile, "tresProfX" );
  addHIST( m_tresProfY, hfile, "tresProfY" );
  addHIST( m_numClsDist, hfile, "numCls" );
  addHIST( m_dirzDist, hfile, "dirZ" );
  addHIST( m_armsDist, hfile, "arms" );
  addHIST( m_lrec, hfile, "lrec" );
  addHIST( m_ldigi, hfile, "ldigi" );
  addHIST( m_lcls, hfile, "lcls" );
  addHIST( m_largeMulGTRC, hfile, "largeMulGTRC" );
  if( iRoot+1 == nRoot ) m_largeMulGTRC->Scale( 1.0/nRoot );
  if( m_badStrips ){
    char hname[]="brms0";
    for( int i=0; i<g_nLayer/3; i++){ 
      sprintf( hname, "brms%d", i );
      addHIST( m_brmsDist[i], hfile, hname );
    }
    addHIST( m_locc, hfile, "locc" );
    addHIST( m_ltrk, hfile, "ltrk" );
    addHIST( m_dist, hfile, "dist" );
    addHIST( m_aPos[0], hfile, "apos0" );
    addHIST( m_aPos[1], hfile, "apos1" );
    addHIST( m_aPos[2], hfile, "apos2" );
    addHIST( m_aPos[3], hfile, "apos3" );
    addHIST( m_aPos[4], hfile, "apos4" );
    //m_occDist->Add( (TH1F*)hfile->FindObjectAny( "occDist" ) );
    //m_poissonDist->Add( (TH1F*)hfile->FindObjectAny( "poissonDist" ) );

    for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw)
      m_towerVar[tw].readHists( hfile, iRoot, nRoot );
  }
  else{
    //m_fracBatTot->Add( (TH1F*)hfile->FindObjectAny( "fracBadTot" ) );
    //m_fracErrDist->Add( (TH1F*)hfile->FindObjectAny( "fracErrDist" ) );
    //m_chisqDist->Add( (TH1F*)hfile->FindObjectAny( "chisqDist" ) );
    //m_chargeScale->Add( (TH1F*)hfile->FindObjectAny( "chargeScale" ) );
    //m_langauWidth->Add( (TH1F*)hfile->FindObjectAny( "langauWidth" ) );
    //m_langauGSigma->Add( (TH1F*)hfile->FindObjectAny( "langauSigma" ) );
    addHIST( m_dirProfile, hfile, "dirProfile" );
    for( int i=4; i!=-1; i--){
      char hname[]="chargeAll";
      if( i<4 ) sprintf( hname, "charge%d", i );
      addHIST( m_chist[i], hfile, hname );
    }
    for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
      int tower = m_towerVar[ tw ].towerId;
      for(int unp = 0; unp != g_nUniPlane; ++unp) {
	layerId lid( unp );
	std::string lname = lid.getLayerName();
	for(int iDiv = 0; iDiv != m_nDiv; ++iDiv){
	  char name[] = "chargeT00X00fe0000";
	  if( m_nDiv != 1 )
	    sprintf(name,"chargeT%d%sfe%d", tower, lname.c_str(), iDiv);
	  else
	    sprintf(name,"chargeT%d%s", tower, lname.c_str());
	  TH1F* hist = (TH1F*)hfile->FindObjectAny( name );
	  for( int ibin=0; ibin!=nTotHistBin; ibin++)
	    m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin] 
	      += (int)hist->GetBinContent( ibin+1 );
	}
      }
    }
  }

  return true;
}


bool totCalib::readBadStripsTxtFiles(const std::string dir, 
				    const std::string names ){  
  std::vector<std::string> fnames;
  splitWords( fnames, names );
  for(std::vector<std::string>::const_iterator fname = fnames.begin();
      fname != fnames.end(); ++fname) {
    if( !readBadStripsTxtFile( dir+'/'+(*fname) ) ) return false;
  }
  return true;
}  

bool totCalib::readBadStripsTxtFile( const std::string filename ){

  if( ! checkFile( filename ) ){
    std::cout << "Invalid text file path: " << filename << std::endl;
    m_log << "Invalid text file path: " << filename << std::endl;
    return false;
  }

  std::cout << "Open bad strips text file: " << filename << std::endl;
  m_log << "Bas strips text file: " << filename << std::endl;

  std::string line, column, word;
  std::vector<std::string> words;
  std::vector<int> towerPtr;
  ifstream is( filename.c_str() );
  getline( is, line );
  UInt_t nl = 0;
  while( is ){
    std::string::size_type pos=0, i=0;
    UInt_t nc = 0;
    int unp = -1;
    while( i != string::npos ){
      i = line.find('\t', pos);
      if(i != string::npos) column = line.substr(pos, i-pos);
      else column = line.substr(pos); // end of line
      pos = i + 1;
      //
      // first line includes module information
      // search for matching module
      //
      if( nl == 0 ){
	std::string serial = "TkrFM" + column;
	int pointer = -1;
	for( UInt_t twr=0; twr<m_towerVar.size(); twr++)
	  if( serial == m_towerVar[twr].hwserial ){ 
	    pointer = twr;
	    std::cout << serial << " " << nc << ";";
	  }
	towerPtr.push_back( pointer );
	continue;
      }
      int twr = -1;
      if( nc < towerPtr.size() ) twr = towerPtr[ nc ];
      if( twr >= 0 ) std::cout << " T" << m_towerVar[twr].towerId;
      //
      // check layer name
      // 
      if( nc == 0 ){
	if( column[0]=='X' or column[0]=='Y' ){
	  std::cout << std::endl << column;
	  int view = 0;
	  if( column[0] == 'Y' ) view = 1;
	  int layer = atoi( column.erase(0,1).c_str() );
	  layerId lid(  layer, view );
	  unp = lid.uniPlane;
	  std::cout  << " u" << unp;
	}
	else unp = -1;
      }
      nc++;
      if( twr < 0 or unp < 0 ) {
	//std::cout << column << ": " << unp << " " << twr << " " << nc << " " << nl;
	//std::cout << std::endl;
	continue;
      }

      //
      // register bad strips
      //
      if( column[0] == '"' ) column.erase(0,1);
      if( column[column.size()-1] == '"' ) column.erase(column.size()-1,1);
      words.clear();
      splitWords( words, column );
      int nst = 0;
      for( UInt_t j=0; j<words.size(); j++){
	word = words[j];
	if( word.size() > 0 ){
	  if( word[word.size()-1] == ',' ) word.erase( word.size()-1,1);
	  int strip = atoi( word.c_str() );
	  m_towerVar[twr].bsVar[unp].knownBadStrips[2].push_back( strip );
	  nst++;
	}
      }
      std::cout << " " << nst;
    }
    nl++;
    getline( is, line );
  }
  is.close();
  
  std::cout << std::endl;
  return true;

}


bool totCalib::readInputXmlFiles(const std::string dir, 
				 const std::vector<std::string>& runIds ){

  for(std::vector<std::string>::const_iterator run = runIds.begin();
      run != runIds.end(); ++run) {
    if( m_badStrips ){
      if( !readBadStripsXmlFile( dir, (*run) ) )
	return false;
    }
    else if( !readTotConvXmlFile( dir, (*run) ) )
      return false;
  }
  return true;
}


bool totCalib::readLatcTfeXmlFile( const std::string filename ){

  if( ! checkFile( filename ) ){
    std::cout << "Invalid latc xml file path: " << filename << std::endl;
    m_log << "Invalid latc xml file path: " << filename << std::endl;
    return false;
  }

  std::cout << "Open LATC TFE xml file: " << filename << std::endl;
  m_log << "LATC TFE xml file: " << filename << std::endl;
  
  using namespace xmlBase;

  XmlParser* parser = new XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse( filename.c_str() );
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    return false;
  }

  if (doc != 0) {  // successful
    //std::cout << "Document successfully parsed" << std::endl;

    // look up generic attributes
    DOMElement* docElt = doc->getDocumentElement();

    // look loop over tower TEM
    DOMNode* temNode = Dom::getFirstChildElement(docElt);
    while( temNode ){ //loop each tem entry.
      int towerId = -1;
      try {
	towerId = Dom::getIntAttribute(temNode, "ID");
      }
      catch (DomException ex) {
	std::cout << "DomException:  " << ex.getMsg() << std::endl;
      }

      if( towerId < 0 ) return false;
      int twr = m_towerPtr[towerId];
      int nChannelFE = g_nStrip / g_nFEC;
      std::string hwserial = m_towerVar[twr].hwserial;
      std::vector<std::string> stripList;
      std::cout << "TEM: " << towerId << ", serial: " << hwserial 
		<< std::endl;

      DOMNode* sptNode = Dom::getFirstChildElement(temNode);
      while( sptNode ){ //loop each layer end entry.
	std::string spt = Dom::getAttribute(sptNode, "ID");
	layerId lid( spt, towerId );
     
	//get first child element
	DOMElement* feNode = Dom::getFirstChildElement(sptNode);
	while( feNode ){
	  int feId = Dom::getIntAttribute(feNode, "ID");
	  DOMElement* maskNode = Dom::findFirstChildByName(feNode,"data_mask");
	  std::string hexmask = Dom::getTextContent(maskNode) ;
	  std::vector<int> chList = bitFlagToList( ~axtoi( hexmask.c_str() ) );
	  //ULong64_t mask = ~axtoi( hexmask.c_str() );
	  //std::cout << towerId << " " << spt << " " << lid.getLayerName() << " " << feId << " " << hexmask << " = " << mask << std::endl;
	  for( UInt_t ich=0; ich!=chList.size(); ich++){
	    int ch = chList[ich];
	    int strip = feId * nChannelFE + ch;
	    //std::cout << feId << " " << ch << " " << strip << std::endl;
	    if( feId < 0 || feId >= g_nFEC || ch < 0 || ch > nChannelFE ){
	      std::cout << towerId << " " << spt
			<< ", invalid strip id: " << strip << std::endl;
	      continue;
	    }
	    // masked strip is considered hot
	    m_towerVar[twr].bsVar[lid.uniPlane].knownBadStrips[0].push_back( strip );
	  }
	  feNode = Dom::getSiblingElement( feNode );
	}
	sptNode = Dom::getSiblingElement( sptNode );
      }
      temNode = Dom::getSiblingElement( temNode );
    }
  }
  else return false;
  delete doc;
  delete parser;
  return true;
}


bool totCalib::readBadStripsXmlFile(const std::string path, 
				    const std::string runid ){
  string filename;
  char fname[] = "/398000364/LICOS/analysis_TEM12/TkrNoiseAndGain_398000364_dead.xml";
  if( atoi(runid.c_str()) < 100000000 ){  // LICOS
    for( int tem=0; tem<g_nTower; tem++){
      sprintf(fname,"/%s/LICOS/analysis_TEM%d/TkrNoiseAndGain_%s_dead.xml", 
	      runid.c_str(), tem, runid.c_str() );
      filename = path + fname;
      //std::cout << filename << std::endl;
      if( ! readBadStripsXmlFile( filename, false ) ) return false;
    }
  }
  else{
    sprintf(fname,"/%s/TkrNoiseAndGain_%s_dead.xml", 
	    runid.c_str(), runid.c_str() );
    filename = path + fname;
    return readBadStripsXmlFile( filename, false );
  }
  return true;
}


bool totCalib::readBadStripsXmlFile(const std::string filename, bool hotStrips ){ 

  if( ! checkFile( filename ) ){
    std::cout << "Invalid bad strips xml file path: " << filename << std::endl;
    m_log << "Invalid bad strips xml file path: " << filename << std::endl;
    return false;
  }

  std::cout << "Open xml file: " << filename << std::endl;
  m_log << "Bad strips xml file: " << filename << std::endl;
  
  //
  // remove dtd from input xml and save to temporary file.
  //
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
  char fname[] = "/tmp/temp-060418-000000.xml";
  sprintf(fname, "/tmp/temp%s.xml", m_timeStamp.c_str() );
  os.open( fname );
  os.write( buffer, length );
  os.close();
  delete buffer;

  using namespace xmlBase;

  XmlParser* parser = new XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse( fname );
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
    return false;
  }
  os.open( fname ); // erase contents of temporary file
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
  delete doc;
  delete parser;
  return true;
}


bool totCalib::readTotConvXmlFile( const std::string path, 
				   const std::string runid ){
  std::string filename;
  char fname[] = "/398000364/LICOS/analysis_TEM12/TkrTotGain_398000364.xml";

  if( atoi(runid.c_str()) < 100000000 ){  // LICOS
    for( int tem=0; tem<g_nTower; tem++){
      sprintf(fname,"/%s/LICOS/analysis_TEM%d/TkrTotGain_%s.xml", 
	      runid.c_str(), tem, runid.c_str() );
      filename = path + fname;
      if( ! checkFile( filename ) ){
	sprintf(fname,"/%s/LICOS/analysis_TEM%d/TkrTotGain_%d.xml", 
		runid.c_str(), tem, atoi(runid.c_str()) );
	filename = path + fname;
      }
      //std::cout << filename << std::endl;
      if( ! readTotConvXmlFile( filename ) ) return false;
    }
  }
  else{
    sprintf(fname,"/%s/TkrTotGain_%s.xml", runid.c_str(), runid.c_str() );
    filename = path + fname;
    //std::cout << filename << std::endl;
    return readTotConvXmlFile( filename );
  }
  return true;
}


bool totCalib::readTotConvXmlFile( const std::string filename ){

  if( ! checkFile( filename ) ){
    std::cout << "Invalid TOT xml file path: " << filename << std::endl;
    m_log << "Invalid TOT xml file path: " << filename << std::endl;
    return false;
  }
  std::cout << "Open xml file: " << filename << std::endl;
  m_log << "TOT xml file: " << filename << std::endl;
  

  using namespace xmlBase;

  XmlParser* parser = new XmlParser(true);
  DOMDocument* doc = 0;
  try {
    doc = parser->parse(filename.c_str());
  }
  catch (ParseException ex) {
    std::cout << "caught exception with message " << std::endl;
    std::cout << ex.getMsg() << std::endl;
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
    if( defElt ){
      std::cout << "INFO: default TOT parameters exist, which are used INSTEAD of those for individual channels." << std::endl;
      m_log << "INFO: default TOT parameters exist, which are used INSTEAD of those for individual channels." << std::endl;
      int tw = m_towerPtr[ towerId ];
      for( int unp=0; unp!=g_nUniPlane; unp++)
	for( int stripId =0; stripId!=g_nStrip; stripId++){
	  m_towerVar[tw].tcVar[unp].totThreshold[stripId] = intercept;
	  m_towerVar[tw].tcVar[unp].totGain[stripId] = slope;
	  m_towerVar[tw].tcVar[unp].totQuad[stripId] = quad;
	}
      delete doc;
      delete parser;
      return true;
    }
    if( len != g_nLayer*g_nView ){
      std::cout << "WARNING: # of layers in xml is not what expected, proceed with care. " << len << std::endl;
      m_log << "WARNING: # of layers in xml is not what expected, proceed with care. " << len << std::endl;
      //return false;
    }

    std::vector<std::string> keywords;
    keywords.push_back( "intercept" );
    keywords.push_back( "slope" );
    keywords.push_back( "quad" );

    int numChannels = 0;
    int numWarn = 0;
    for(int i=0; i<len; i++){ // loop for each uniplane entry
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
      while( getParam( elder, lid, keywords ) ){
	younger = Dom::getSiblingElement( elder );
	elder = younger;
	numStrip++;
	numChannels++;
      }
      if( numStrip != g_nStrip && numWarn < 5 ){
	std::cout << "WARN: # of strips in uniplane is not what expected, " 
		  << numStrip << std::endl;
	m_log << "WARN: # of strips in uniplane is not what expected, " 
	      << numStrip << std::endl;
	numWarn++;
	if( numWarn == 5 ){
	  std::cout << "WARN: this is the last warning, further warnings are suppressed." << std::endl;
	  m_log << "WARN: this is the last warning, further warnings are suppressed." << std::endl;
	}
	//return false;
      }
    }
    if( numChannels != g_nLayer*g_nView*g_nStrip ){
	std::cout << "ERROR: # of strips in xml is invalid, " 
		  << numChannels << std::endl;
	m_log << "ERROR: # of strips in xml is invalid, " 
	      << numChannels << std::endl;
	return false;
    }

  }
  else return false;

  delete parser;
  delete doc;
  return true;
}

bool totCalib::getParam(const DOMElement* totElement, layerId lid, std::vector<std::string> keywords ){  

  using namespace xmlBase;

  int tw = m_towerPtr[ lid.tower ];
  int unp = lid.uniPlane;
  std::string lname = lid.getLayerName();

  int stripId;
  std::vector<float> values;

  try{
    stripId = Dom::getIntAttribute( totElement, "id" ); 
  } //if there isn't next strip,go to next layer or view
  catch(DomException ex){
    //cout << "finished " << lname << endl;
    return false;
  }
  if( stripId < 0 || stripId >= g_nStrip ){
    std::cout << "ERROR: " << lname
	      << ", Invalid strip id: " << stripId << std::endl;
    m_log << "ERROR: " << lname 
	  << ", Invalid strip id: " << stripId << std::endl;
    return true;
  }

  for( unsigned int i=0; i<keywords.size(); i++){
    try{
      float value = Dom::getDoubleAttribute( totElement, keywords[i].c_str() );
      values.push_back( value );
    }
    catch(DomException ex){
      cout << "ERROR, no attribute for " << keywords[i] << ": (L,S)=(" 
	   << lname << ", "  << stripId << ")" << endl;
      m_log << "ERROR, no attribute for " << keywords[i] << ": (L,S)=(" 
	    << lname << ", "  << stripId << ")" << endl;
    }
  }
  /*  cout <<"stripId" << stripId
      <<",offset" << offset
      <<",gain" << gain
      <<",quad" << quad <<endl;*/
  
  m_towerVar[tw].tcVar[unp].totThreshold[stripId] = values[0];
  m_towerVar[tw].tcVar[unp].totGain[stripId] = values[1];
  m_towerVar[tw].tcVar[unp].totQuad[stripId] = values[2];
  return true;
}


float totCalib::calcCharge( layerId lid, int iStrip, int tot) const
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  int tw = m_towerPtr[ lid.tower ];
  int unp = lid.uniPlane;
  // TOT to charge conversion
  float charge = m_towerVar[tw].tcVar[unp].totThreshold[iStrip] 
    + time*m_towerVar[tw].tcVar[unp].totGain[iStrip]
    + time*time*m_towerVar[tw].tcVar[unp].totQuad[iStrip];
  
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
  
  for( int unp = 0; unp < g_nUniPlane; unp++){
    layerId lid( unp );
    int tray = lid.tray;
    std::string which = lid.which;
    std::string lname = lid.getLayerName();

    xmlFile << std::endl
	    << "   <!-- **** layer " << lname  
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
  

void totCalib::findBadStrips()
{

  float sumThreshold, occThreshold;

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
      int unp = uniPlane;
      std::string lname = lid.getLayerName();

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
      std::vector<float> vratio;
      float ratio;
      for( int strip=0; strip!=g_nStrip; strip++){
	if( m_towerVar[tw].bsVar[uniPlane].lHits[strip] == 0 ) continue;
	if( m_towerVar[tw].bsVar[uniPlane].tHits[strip] == 0 ) continue;
	ratio = (float)m_towerVar[tw].bsVar[uniPlane].lHits[strip] 
	  / (m_towerVar[tw].bsVar[uniPlane].tHits[strip]) ;
	vratio.push_back( ratio );
      }
      float eff = getTrunctedMean( vratio, 0.9 );
      m_log << "T" << tower << " " << lname << ": " << eff;
      if( eff < minEff ) m_log << "*****";
      m_log << std::endl;
	
      //
      // main bad strips finder
      //
      int numDeadStrips = 0;
      float meanRatio = 0.0;
      int numGood = 0;
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
	  if(lid.layer%2 == 0) jWafer = g_nWafer-1-iWafer;
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

      m_log << "T" << tower << " " << lname << ", # of bad channels:";
      std::cout << "T" << tower << " " << lname << ", # of bad channels:";
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
  std::string lname = lid.getLayerName();
  //
  // initialize known bad strip flag
  //
  bool known[g_nStrip], dead[g_nStrip], hot[g_nStrip], conflict[g_nStrip],
    found[g_nStrip], debug=false;
  for( int i=0; i!=g_nStrip; i++){
    known[i] = false;
    dead[i] = false;
    hot[i] = false;
    conflict[i] = false;
    found[i] = false;
  }
  //
  // register bad strips found in this job (iBad=g_nBad-1)
  //
  int ist;
  if( m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].size() != 0 ){
    for( UInt_t i=0; i!=m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].size(); i++){
      ist = m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1][i];
      if( found[ist] ){
	std::cout << "Redundant bad strip: T" << lid.tower << " " 
		  << lname << " " << ist << std::endl;
	m_log << "Redundant bad strip: T" << lid.tower << " " 
	      << lname << " " << ist << std::endl;
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
	// disappeared disocnnected strip, categorize as intermittent
	m_towerVar[tw].bsVar[unp].badStrips[4].push_back( ist );
	m_towerVar[tw].bsVar[unp].badStrips[g_nBad-1].push_back( ist );
	std::cout << "Missing known bad strip: T" << lid.tower << " " 
		  << lname << " " << ist << std::endl;
	m_log << "Missing known bad strip: T" << lid.tower << " " 
	      << lname << " " << ist << std::endl;
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
		  << lname << " " << ist << std::endl;
	m_log << "Missing known dead strip: T" << lid.tower << " " 
	      << lname << " " << ist << std::endl;
	debug = true;
      }
      if( iBad == 0 ){
	hot[ist] = true;
	if( known[ist] ){
	  std::cout << "Conflict strip: T" << lid.tower << " " 
		    << lname << " " << ist << std::endl;
	  m_log << "Conflict strip: T" << lid.tower << " " 
		<< lname << " " << ist << std::endl;
	  known[ist] = false;
	  conflict[ist] = true;
	  m_towerVar[tw].bsVar[unp].badStrips[6].push_back( ist );
	}
      }
      else dead[ist] = true;
      m_towerVar[tw].bsVar[unp].badStrips[iBad].push_back( ist );
    }
    //
    // print debug information when inconsistency is found.
    //
    if( debug ){
      std::cout << "DEBUG: T" << lid.tower << " " << lname
		<< ", " << iBad << ": ";
      m_log << "DEBUG: T" << lid.tower << " " << lname
	    << ", " << iBad << ": ";
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
      //std::cout << unp << " " << ist << ", " << known[ist] << found[ist] << dead[ist] << hot[ist] << std::endl;
      //m_log << unp << " " << ist << ", " << known[ist] << found[ist] << dead[ist] << hot[ist] << std::endl;
      if( known[ist] ) continue; // leave disconnected strips found online
      if( !found[ist] ){
	std::cout << "Inconsistent dead strip: T" << lid.tower << " " 
		  << lname << " " << ist << std::endl;
	m_log << "Inconsistent bad strip: T" << lid.tower << " " 
	      << lname << " " << ist << std::endl;
	m_towerVar[tw].bsVar[unp].knownBadStrips[g_nBad-1].push_back( ist );
	debug = true;
      }
      if( (iBad<g_nBad-2 && (dead[ist] || hot[ist])) ){ 
	m_towerVar[tw].bsVar[unp].badStrips[iBad][i] = g_nStrip;
	numMatch++;
      }
    }
    //
    // print debug information when inconsistency is found.
    //
    if( debug ){
      std::cout << "DEBUG: T" << lid.tower << " " << lname 
		<< ", " << iBad << ": ";
      m_log << "DEBUG: T" << lid.tower << " " << lname
	    << ", " << iBad << ": ";
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
		    << lname << " " << ist << std::endl;
	  m_log << "Inconsistent matched strip: T" << lid.tower << " " 
		<< lname << " " << ist << std::endl;
	  std::cout << "DEBUG: T" << lid.tower << " " << lname
		    << ", " << iBad << " " << numMatch << ": ";
	  m_log << "DEBUG: T" << lid.tower << " " << lname
		<< ", " << iBad << " " << numMatch << ": ";
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

  std::ofstream latxml, lat2xml, fmxml;
  if( m_towerVar.size() > 1 ){
    sprintf( fname, "/LAT_DeadStrips_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    filename += fname;

    latxml.open( filename.c_str() );
    if( latxml ){
      std::cout << "Open LAT dead strips xml file: " << filename << std::endl;
      m_log << "LAT dead strips xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }
    openBadStripsXml( latxml, dtd );

    //
    // file for bad strips xml file which includes hot strips as well.
    //
    filename = m_outputDir;
    sprintf( fname, "/LAT_BadStrips_%s-%s.xml", 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );
    filename += fname;

    lat2xml.open( filename.c_str() );
    if( lat2xml ){
      std::cout << "Open LAT bad strips xml file: " << filename << std::endl;
      m_log << "LAT bad strips xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }
    openBadStripsXml( lat2xml, dtd );

  }

  for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
    filename = m_outputDir;
    sprintf( fname, "/%s_DeadStrips_%s-%s.xml", 
	     m_towerVar[tw].hwserial.c_str(), 
	     m_dateStamp.c_str(), m_timeStamp.c_str() );  
    filename += fname;

    fmxml.open( filename.c_str() );
    if( fmxml ){
      std::cout << "Open dead strips xml file: " << filename << std::endl;
      m_log << "Dead strips xml file: " << filename << std::endl;
    }
    else{
      std::cout << filename << " cannot be opened." << std::endl;
      return;
    }

    openBadStripsXml( fmxml, dtd );
    fillTowerBadStrips( fmxml, tw, false, g_nBad ); // add all bad sttrips field
    if( m_towerVar.size() > 1 ){
      fillTowerBadStrips( latxml, tw );
      fillTowerBadStrips( lat2xml, tw, true );
    }

    fmxml << "</badStrips>" << std::endl;
    fmxml.close();
  }

  if( m_towerVar.size() > 1 ){
    latxml << "</badStrips>" << std::endl;
    latxml.close();
    lat2xml << "</badStrips>" << std::endl;
    lat2xml.close();
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
  
  char runs[] = "012345678-012345678";
  sprintf( runs, "%09d-%09d", m_first_run, m_last_run );
  xmlFile << "  <generic calType=\"stripOccupancy\" creatorName=\"badStrips\""
	  << " creatorVersion =\"" << m_version
	  << "\" fmtVersion=\"NA\" instrument=\"TWR\" runId=\"" << runs
	  << "\" timestamp=\"" << m_dateStamp << m_timeStamp << "\">" << std::endl
	  << "    <inputSample mode=\"NA\" source=\"CosmicMuon\" startTime=\"" 
	  << m_startTime << "\" stopTime=\"" << m_stopTime 
	  << "\" triggers=\"TKR\">" << std::endl
	  << " Cosmic ray muon data for occupancy analysis " << std::endl
	  << "    </inputSample>" << std::endl
	  << "  </generic>" << std::endl;
  
}
 

void totCalib::fillTowerBadStrips( std::ofstream &xmlFile, const int tw, 
				   const bool saveHotStrips,
				   const int nBad ){
  
  //dead, disconnected, partial disconnected, intermittently disconnected, 
  //intermittently partial disconnexcted
  int howBad[g_nBad] = {1,2,4,12,20,28,64,128}; 
  std::string cBad[g_nBad] = {"hot","dead","disconnected",
			      "partially disconnected",
			      "intermittently disconnected",
			      "intermittently partially connected",
			      "conflict",
			      "all bad (hot/online only included)"};
  
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
    int tray = lid.tray;
    std::string which = lid.which;
    std::string lname = lid.getLayerName();
    
    xmlFile << std::endl
	    << "    <!-- layer " << lname << " -->" << std::endl;
    
    for( int iBad=0; iBad!=nBad; iBad++ ){ 
      int itr = m_towerVar[tw].bsVar[uniPlane].badStrips[iBad].size();
      xmlFile << "    <!-- # of " << cBad[iBad] << " strips: " << itr 
	      << " -->" << std::endl;
      if( !saveHotStrips && iBad == 0 ) 
	xmlFile << "    <!-- " << std::endl; // hot strip is not included.
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
      if( !saveHotStrips && iBad == 0 ) 
	xmlFile << "    --> " << std::endl; // hot strip is not included.
      
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
    lfac += log(float(k));
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
    return log10( prob ) - mean*log10( exp(1.) );
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

