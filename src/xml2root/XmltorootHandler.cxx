//   $Header: /nfs/slac/g/glast/ground/cvs/calibGenTKR/src/xml2root/XmltorootHandler.cxx,v 1.1 2005/03/31 23:13:51 jrb Exp $
/**
   @file XmltorootHandler.cxx

   Implementation of callbacks for Xerces SAX parser
*/

#include <iostream>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "XmltorootHandler.h"
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/AttributeList.hpp>

#include "commonRootData/idents/TkrId.h"
#include "calibRootData/Tkr/TkrTower.h"
#include "calibRootData/Tkr/Tot.h"
#include "calibRootData/Tkr/ChargeScale.h"
#include "facilities/Util.h"

namespace {
  static XMLCh* totX = 0;
  static XMLCh* thresholdX = 0;
  static XMLCh* chargeScaleX = 0;
  static XMLCh* genericX = 0;
  static XMLCh* towerX = 0;
  static XMLCh* uniplaneX = 0;
  static XMLCh* stripTotX = 0;
  static XMLCh* stripScaleX = 0;
  static XMLCh* stripTrgThreshX = 0;
  static XMLCh* stripDataThreshX = 0;
  static XMLCh* gtfeScaleX = 0;
  static XMLCh* gtfeTrgThreshX = 0;
  static XMLCh* gtfeDataThreshX = 0;
  static XMLCh* inputSampleX = 0;
  void initElementNames() {
    totX = XMLString::transcode("tot");
    thresholdX = XMLString::transcode("threshhold");
    chargeScaleX = XMLString::transcode("chargeScale");
    genericX = XMLString::transcode("generic");
    towerX = XMLString::transcode("tower");
    uniplaneX = XMLString::transcode("uniplane");
    stripTotX = XMLString::transcode("stripTot");
    stripScaleX = XMLString::transcode("stripScale");
    stripTrgThreshX = XMLString::transcode("stripTrgThresh"); 
    stripDataThreshX = XMLString::transcode("stripDataThresh"); 
    gtfeScaleX = XMLString::transcode("gtfeScale");
    gtfeTrgThreshX = XMLString::transcode("gtfeTrgThresh");
    gtfeDataThreshX = XMLString::transcode("gtfeDataThresh");
    inputSampleX = XMLString::transcode("inputSample");
  }
}

XmltorootHandler::XmltorootHandler(std::string outname)
  : m_tfile(0), m_treeTower(0), m_totUni(0), m_scaleUni(0), m_uniBranch(0), 
    m_maxStripId(0), m_maxGtfeId(0), m_stripCount(0), m_gtfeCount(0)
{
  m_outname = outname;

  m_tfile = new TFile(outname.c_str(), "RECREATE");
  if (!m_tfile->IsOpen() ) {
    std::cout << "Could not open file " << outname
              << " for writing " << std::endl;
  }
  else {
    std::cout << "Successfully opened " << outname << " for write"
              << std::endl;
  }
  m_tfile->cd();
  if (totX == 0) initElementNames();
}

XmltorootHandler::~XmltorootHandler() {
  // Delete anything needing deletion
}

void XmltorootHandler::endDocument() {
  int nBytes =  m_tfile->Write(0, TObject::kOverwrite);
  std::cout << "Wrote " << nBytes << " to file " << m_outname << std::endl;

  m_tfile->Close();
  delete m_tfile;
  m_tfile = 0;
}

void XmltorootHandler::endElement(const XMLCh* const name) {

  if (XMLString::equals(name, totX)) {
    // maybe need to do something
    return;
  }  else if (XMLString::equals(name, thresholdX)) {
    // will be same as totX, if there is anything at all to do
    return;

  }  else if (XMLString::equals(name, chargeScaleX)) {
    // will be same as totX, if there is anything at all to do
    return;

  }  else if (XMLString::equals(name, towerX)) {
    if (m_totUni) {
      delete m_totUni;          // assuming cal_type is "tot"
      m_totUni = 0;
    }
    if (m_scaleUni) {
      delete m_scaleUni;
      m_scaleUni = 0;
    }
    m_maxStripId = 0;
    m_maxGtfeId = 0;

  }  else if (XMLString::equals(name, uniplaneX)) {
    m_uniBranch->Fill();
    if (m_caltype == "tot") {
      m_totUni->Clear();  
    }
    else if (m_caltype == "scale") {
      m_scaleUni->Clear();
    }

  }  else {  // probably nothing to do for any of the following:
    // generic, cuts, input sample, any of individual data elements for
    // strips or gtfe's
    return;
  }
}

void XmltorootHandler::characters(const XMLCh* const /* chars */, 
                                  const unsigned int /* length */) {
  // don't do anything.  We don't have text content anywhere that matters.
}
void XmltorootHandler::ignorableWhitespace(const XMLCh* const /* chars */,
                                           const unsigned int /* length */) {
  // don't need to do anything
}

void XmltorootHandler::processingInstruction(const XMLCh* const /* target */,
                                             const XMLCh* const  /*  data */) {
  // nor here
}

void XmltorootHandler::startDocument(){
  // probably do nothing here as well

}
              
void XmltorootHandler::startElement(const XMLCh* const name, 
                                   AttributeList& attributes) {
  if (XMLString::equals(name, totX)) {
    // no attributes; just need to keep track of file type
    m_caltype = std::string("tot");
    return;
  }  
  else if (XMLString::equals(name, thresholdX)) {
    m_caltype = std::string("threshold");   // maybe should throw NYI instead
    return;
    
  }  
  else if (XMLString::equals(name, chargeScaleX)) {
    m_caltype = std::string("scale"); 
    return;

  }  
  else if (XMLString::equals(name, genericX)) {
    return;   // for now, anyway
  }  
  else if (XMLString::equals(name, towerX)) {
    return startTower(attributes);
  }  
  else if (XMLString::equals(name, uniplaneX)) {
    return startUniplane(attributes);
  }  
  else if (XMLString::equals(name, stripTotX)) {
    return startStripTot(attributes);

  }  
  else if (XMLString::equals(name, stripScaleX)) {
    return startScaleObj(attributes, true);
  }  
  else if (XMLString::equals(name, stripTrgThreshX)) {
    std::cerr << "Thresholds calibrations not supported.  Exiting.."
              << std::endl;
    std::cerr.flush();
    exit(0);
  }  
  else if (XMLString::equals(name, stripDataThreshX)) {
    std::cerr << "Thresholds calibrations not supported.  Exiting.."
              << std::endl;
    std::cerr.flush();
    exit(0);
  }  
  else if (XMLString::equals(name, gtfeScaleX)) {
    return startScaleObj(attributes, false);
  }  
  else if (XMLString::equals(name, inputSampleX)) {
    return;
  }  
  else if (XMLString::equals(name, gtfeTrgThreshX)) {
    std::cerr << "Thresholds calibrations not supported.  Exiting.."
              << std::endl;
    std::cerr.flush();
    exit(0);
  }   
  else if (XMLString::equals(name, gtfeDataThreshX)) {
    std::cerr << "Thresholds calibrations not supported.  Exiting.."
              << std::endl;
    std::cerr.flush();
    exit(0);
  }   
  else { // a tag we don't care about, at least for now
    // includes cuts, inputSample
    char* tagname = XMLString::transcode(name);

    std::cout << "Encountered unhandled element " << tagname << std::endl;
    XMLString::release(&tagname);
    return;
  }
}

// <tower>       start a tree.  Attributes give us row, col and hwserial
//               so can also make tower branch.  maxStripId 
//               is used in making uniplane branch
void XmltorootHandler::startTower(AttributeList& attributes) {
  // Fetch all the attributes we need to make tree, tower branch
  const XMLCh* rowX = attributes.getValue("row");
  const XMLCh* colX = attributes.getValue("col");
  const XMLCh* serX = attributes.getValue("hwserial");
  //... and unilayer branch.
  const XMLCh* maxStripIdX = attributes.getValue("maxStripId");
  const XMLCh* maxGtfeIdX = attributes.getValue("maxGtfeId");

  char* trans =  XMLString::transcode(rowX);
  std::string rowS = std::string(trans);
  int towerRow = facilities::Util::stringToInt(rowS);
  m_towerRow = towerRow;
  // Make a string that looks like "TowerAB" where A = row#, B = col#
  std::string towS = std::string("Tower") + rowS;
  std::string treeLabel = m_caltype + std::string("data for ") + towS;

  XMLString::release(&trans);
  trans = XMLString::transcode(colX);
  std::string colS = std::string(trans);
  int towerCol  = facilities::Util::stringToInt(colS);
  m_towerCol = towerCol;
  towS += colS;

  XMLString::release(&trans);

  // Even if this is not the first tree to be made, we don't have
  // to worry about saving the pointer to the old one.  We're all
  // done with specific references to it.  Deleting the TFile it
  // belongs to will delete it.
  m_treeTower = new TTree(towS.c_str(), treeLabel.c_str());

  // Uncommenting the following apparently causes garbage to be
  // written to the file rather than the contents of the tree
  //  m_treeTower->SetDirectory(m_tfile);

  std::string serS("");
  if (serX) {
    trans = XMLString::transcode(serX);
    serS = std::string(trans);
    XMLString::release(&trans);
  }
  // Make and fill the tower branch within the tree
  calibRootData::TkrTower tower((unsigned) towerRow, (unsigned) towerCol,
                                TString(serS.c_str()) );
  calibRootData::TkrTower* pTower = &tower;
  char* className = "calibRootData::TkrTower";
  TBranch* towerBranch = m_treeTower->Branch(className, className, &pTower);
  towerBranch->Fill();

  // Also make branch for unilayer objects
  int maxStripId = 1535;  // only need this if file has deprecated dtd
  int maxGtfeId = 23;

  if (maxStripIdX) {
    trans = XMLString::transcode(maxStripIdX);
    std::string maxS = std::string(trans);
    XMLString::release(&trans);
    maxStripId = facilities::Util::stringToInt(maxS);
  }
  m_maxStripId = (unsigned) maxStripId;

  if (maxGtfeIdX) {
    trans = XMLString::transcode(maxGtfeIdX);
    std::string maxS = std::string(trans);
    XMLString::release(&trans);
    maxGtfeId = facilities::Util::stringToInt(maxS);
  }
  m_maxGtfeId = (unsigned) maxGtfeId;


  commonRootData::TkrId nullId;

  if (m_caltype == std::string("tot")) {

    // Unilayer object will have to stick around; can't allocate on stack
    m_totUni = new calibRootData::TotUnilayer(nullId, m_maxStripId + 1);

    m_uniBranch = m_treeTower->Branch("calibRootData::TotUnilayer",
                                      "calibRootData::TotUnilayer",
                                      &m_totUni);
  } else if (m_caltype == std::string("scale")) {
    m_scaleUni = 
      new calibRootData::ChargeScaleUnilayer(nullId);
    m_uniBranch = m_treeTower->Branch("calibRootData::ChargeScaleUnilayer",
                                      "calibRootData::ChargeScaleUnilayer",
                                      &m_scaleUni);
  }
  // otherwise NYI
}


// <uniplane>  Get tray, bot/top attributes; generate TkrId;
//
void XmltorootHandler::startUniplane(AttributeList& attributes) {
  const XMLCh* trayX = attributes.getValue("tray");
  const XMLCh* topX = attributes.getValue("which");

  bool top;
  // set them from trayX, topX above

  m_replicate = false;                         // normal case
  m_stripCount = 0; 
  m_gtfeCount = 0; 

  char* trans = XMLString::transcode(trayX);
  std::string trayS = std::string(trans);
  int tray = facilities::Util::stringToInt(trayS);
  XMLString::release(&trans);

  trans = XMLString::transcode(topX);
  top = (strcmp(trans, "top") == 0);
  XMLString::release(&trans);
  
  // Constructor for TkrId has arguments towerX, towerY, corresponding
  // to column, row
  if (m_caltype == "tot") {
    m_totUni->setId(commonRootData::TkrId(m_towerCol, m_towerRow, tray, top));
  }
  else if (m_caltype == "scale") {
    m_scaleUni->setId(commonRootData::TkrId(m_towerCol, m_towerRow, tray, 
                                            top));
  }  
}

// <stripTot>  
void XmltorootHandler::startStripTot(AttributeList& attributes) {
  if (m_replicate) {               // this is an error
    std::cerr << "Cannot handle additional strips when m_replicate is 'true' " << std::endl;
    return;
  }

  const XMLCh* idX = attributes.getValue("id");
  const XMLCh* slopeX = attributes.getValue("slope");
  const XMLCh* interceptX = attributes.getValue("intercept");
  const XMLCh* quadX = attributes.getValue("quad");
  const XMLCh* chi2X = attributes.getValue("chi2");
  const XMLCh* dfX = attributes.getValue("df");

  if ((!idX) || (!slopeX) || (!interceptX) || (!quadX) || (!chi2X) 
      || (!dfX) ) {
    std::cerr << "Missing attribute values for <stripTot>" << std::endl;
    std::cerr << "Are you using the correct dtd?    Exiting..." << std::endl;
    std::cerr.flush();
    exit(0);
  }

  char* trans = XMLString::transcode(idX);
  std::string s = std::string(trans);
  int id = facilities::Util::stringToInt(s);
  XMLString::release(&trans);

  if ((unsigned) id > 1+ m_maxStripId) {
    std::cerr << "Cannot include out-of range strip id " << id 
             << " where max is  " << m_maxStripId << std::endl;
    return;
  } 
  if ((unsigned) id == (m_maxStripId + 1) ) {
    if (m_stripCount > 0) {  // error
      std::cerr << "Can't have more than one strip when replicating" 
                << std::endl;
      return;
    }
    m_replicate = true;
  }
  else m_stripCount++;   // normal case; handle 1 strip

  trans = XMLString::transcode(slopeX);
  s = std::string(trans);
  float slope = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(interceptX);
  s = std::string(trans);
  float intercept = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(quadX);
  s = std::string(trans);
  float quad = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(chi2X);
  s = std::string(trans);
  float chi2 = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(dfX);
  s = std::string(trans);
  float df = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  if (!m_replicate) {
    calibRootData::TotStrip strip(id, slope, intercept, quad, chi2, df);
    m_totUni->setStrip(strip);
  }
  else {
    for (unsigned i = 0; i <= m_maxStripId; i++) {
      calibRootData::TotStrip strip(i, slope, intercept, quad, chi2, df);
      m_totUni->setStrip(strip);
    }
  }
}
void XmltorootHandler::startScaleObj(AttributeList& attributes, bool isStrip) {
  if (m_replicate) {               // this is an error
    std::cerr << 
    "Cannot handle additional child scale objects when m_replicate is 'true' " 
              << std::endl;
    return;
  }
  unsigned* toIncr = 0;
  unsigned maxId;
  if (isStrip) {
    toIncr = &m_stripCount;
    maxId = m_maxStripId;
  }
  else {
    toIncr = &m_gtfeCount;
    maxId = m_maxGtfeId;
  }

  if (*toIncr == 0) { // finish setting up m_scaleUni
    m_scaleUni->resize(maxId+1, isStrip);
  }
  const XMLCh* idX = attributes.getValue("id");
  const XMLCh* scaleX = attributes.getValue("chargeScale");
  const XMLCh* errorX = attributes.getValue("error");
  const XMLCh* chi2X = attributes.getValue("chi2");
  const XMLCh* dfX = attributes.getValue("df");
  
  if ((!idX) ||  (!scaleX) || (!errorX) || (!chi2X) 
      || (!dfX) ) {
    std::cerr << "Missing attribute values for <stripScale> or <gtfeScale>" 
              << std::endl;
    std::cerr << "Are you using the correct dtd?    Exiting..." << std::endl;
    std::cerr.flush();
    exit(0);
  }


  char* trans = XMLString::transcode(idX);
  std::string s = std::string(trans);
  int id = facilities::Util::stringToInt(s);
  XMLString::release(&trans);
  
  if ((unsigned) id > 1+ maxId) {
    std::cerr << "Cannot include out-of range  id " << id 
              << " where max is  " << maxId << std::endl;
    return;
  } 
  // look for special case, telling us to replicate for default calib.
  if ((unsigned) id == (maxId + 1) ) { 
    if (*toIncr > 0) {  // error
      std::cerr << "Can't have more than one strip when replicating" 
                << std::endl;
      return;

    }
    m_replicate = true;
  }
  else (*toIncr)++;   // normal case; handle 1 strip
  
  trans = XMLString::transcode(scaleX);
  s = std::string(trans);
  float scale = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);
  
  trans = XMLString::transcode(errorX);
  s = std::string(trans);
  float error = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(chi2X);
  s = std::string(trans);
  float chi2 = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  trans = XMLString::transcode(dfX);
  s = std::string(trans);
  float df = facilities::Util::stringToDouble(s);
  XMLString::release(&trans);

  if (!m_replicate) {
    calibRootData::ChargeScaleObj 
      child(id, scale, error, chi2, df);
    m_scaleUni->setChild(child);
  }
  else {
    for (unsigned i = 0; i <= maxId; i++) {
      calibRootData::ChargeScaleObj child(i, scale, error, chi2, df);
      m_scaleUni->setChild(child);    
    }
  }
}   // end of startStripScale 
