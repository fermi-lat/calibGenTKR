#ifndef XmltoorootHandler_h
#define XmltoorootHandler_h
// $Header: /nfs/slac/g/glast/ground/cvs/users/jrb/xmlToRoot/src/XmltorootHandler.h,v 1.2 2005/03/28 20:01:02 jrb Exp $
#include    <xercesc/sax/HandlerBase.hpp>
#include    <xercesc/sax/AttributeList.hpp>
#include    <string>

XERCES_CPP_NAMESPACE_USE

class TFile;
class TTree;
class TBranch;

namespace calibRootData {
  class TotUnilayer;
  class ChargeScaleUnilayer;
}

/**
   @class XmltorootHandler

   Defines callbacks needed by SAX DocumentHandler interface, as well as 
   private data  and routines needed for converting TKR XML calibration files
   to ROOT.  Currently handles ToT gain and Charge scale.

   At this time the only non-null implementations among the DocumentHandler 
   interface routines are startElement, endElement, and endDocument.
*/
class XmltorootHandler : public HandlerBase {
public:
  XmltorootHandler(std::string outname);
  ~XmltorootHandler();

  /// Implementation needed for SAX DocumentHandler interface
  void endDocument();

  /// Implementation needed for SAX DocumentHandler interface
  void endElement(const XMLCh* const name);

  /// Implementation needed for SAX DocumentHandler interface
  void characters(const XMLCh* const chars, const unsigned int length);

  /// Implementation needed for SAX DocumentHandler interface
  void ignorableWhitespace
  (
   const   XMLCh* const    chars
   , const unsigned int    length
   );

  /// Implementation needed for SAX DocumentHandler interface
  void processingInstruction
  (
   const   XMLCh* const    target
   , const XMLCh* const    data
   );

  /// Implementation needed for SAX DocumentHandler interface
  void startDocument();

  /// Implementation needed for SAX DocumentHandler interface
  void startElement(const XMLCh* const name, AttributeList& attributes);


private:
  std::string m_outname;

  /// Current TFile  
  TFile*  m_tfile;

  /// Tree for tower in progress
  TTree*  m_treeTower;

  /// Keep information on current uniplane for ToT file
  calibRootData::TotUnilayer* m_totUni;

  /// Keep information on current uniplane for charge scale file
  calibRootData::ChargeScaleUnilayer* m_scaleUni;
  TBranch* m_uniBranch;

  /// one of "tot", "threshold", "scale"
  std::string m_caltype;  

  // Keep around for generating commoRootData::TkrId
  unsigned m_towerRow;    
  unsigned m_towerCol;

  /// flag set to true upon recognizing signal in xml input; 
  /// cleared for each uni.
  bool m_replicate;  

  // Keep around so that we know how big a TClonesArray to allocate,
  // if we need one.
  unsigned m_maxStripId;
  unsigned m_maxGtfeId;
  unsigned m_stripCount;  // how many have we seen so far?
  unsigned m_gtfeCount;  // how many have we seen so far?

  /// Called from startElement to handle <tower>
  void startTower(AttributeList& attributes);
  /// Called from startElement to handle <uniplane>
  void startUniplane(AttributeList& attributes);
  /// Called from startElement to handle <stripTot>
  void startStripTot(AttributeList& attributes);
  /// Called from startElement to handle <stripScale> or <stripGtfe>
  void startScaleObj(AttributeList& attributes, bool isStrip);

  // Will need to keep track of active tree, active branch,...  
};

#endif
