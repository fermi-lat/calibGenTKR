/**
  $Header: /nfs/slac/g/glast/ground/cvs/users/jrb/xmlToRoot/src/xmltoroot.cxx,v 1.2 2005/03/28 20:01:02 jrb Exp $
 stand-alone reads in a Tracker XML calibration file (e.g. ToT signal)
 and writes out equivalent ROOT file. 

 Use Xerces sample program SAXPrint as partial model.
*/

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include "xmlBase/XmlErrorHandler.h"

//// #include "SAXPrint.hpp"
// instead will need something like..
#include <string>
#include <iostream>

// ..and a #include for handler .h file
#include "XmltorootHandler.h"

// For now don't let caller change any of these things.  These
// settings are probably fine for our limited application.
static bool                     doNamespaces        = false;
static bool                     doSchema            = false;
static bool                     schemaFullChecking  = false;
// static const char*              encodingName    = "LATIN1";
static char*                    xmlFile         = 0;
static SAXParser::ValSchemes    valScheme       = SAXParser::Val_Auto;


/// Print out usage information
static void usage()
{
  //  XERCES_STD_QUALIFIER cout 
  std::cout 
    << "\nUsage:\n"
    "    xmltoroot <XML input file>   [<ROOT output file>] \n\n"
    "This program reads in an xml calibration file of certain\n"
    "supported types and writes out a corresponding ROOT file\n"
    "ROOT output filename is optional. If omitted, defaults to\n"
    "input filename with extension = .root\n"
    <<  std::endl;
    //    <<  XERCES_STD_QUALIFIER endl;
}


int main(int argC, char* argV[])
{
  // Initialize Xerces
  try
  {
    XMLPlatformUtils::Initialize();
  }

  catch (const XMLException& toCatch)
  {
    char* charMsg = XMLString::transcode(toCatch.getMessage());
    std::cerr << "Error during initialization! :\n"
              << charMsg
              << std::endl;
    XMLString::release(&charMsg);
    return 1;
  }

  // Check command line and extract arguments.
  if (argC < 2)
  {
    usage();
    XMLPlatformUtils::Terminate();
    return 1;
  }

  // Watch for special case help request
  if ( (!strcmp(argV[1], "-?")) ||
       (!strcmp(argV[1], "-help")) ||
       (!strcmp(argV[1], "-h")) ||
       (!strcmp(argV[1], "--help")) )
  {
    usage();
    XMLPlatformUtils::Terminate();
    return 2;
  }

  // otherwise translate a file
  std::string rootFile("");

  xmlFile = argV[1];
  int errorCount = 0;

  if (argC > 2) {
    rootFile = std::string(argV[2]);
  }
  else {
    // xml file spec had better end with .xml
    rootFile = std::string(xmlFile);
    unsigned xnameSize = rootFile.size();
    std::string ftype;
    //    ftype.assign(rootFile, 0, xnameSize-4);
    ftype.assign(rootFile, xnameSize-4, 4);
    if (ftype != std::string(".xml")) {  // error
      std::cout << "Specify 2nd arg. for output file or use .xml extension\n";
      std::cout << "on input file." << std::endl;
      return 3;
    }
    else {
      //      rootFile = rootFile.subString(0, xnameSize - 3) + "root";
      rootFile.resize(xnameSize - 3);
      rootFile += "root";
    }

  }
  //
  //  Create a SAX parser object. For the time being, maybe forever,
  //  don't give user any way to change these settings.
  //  Might want to make wrapper class in xmlBase for SAXParser
  //  as was done for DOMParser.
  //
  SAXParser* parser = new SAXParser;
  parser->setValidationScheme(valScheme);
  parser->setDoNamespaces(doNamespaces);
  parser->setDoSchema(doSchema);
  parser->setValidationSchemaFullChecking(schemaFullChecking);


  //   To be dealt with:


  //
  //  xmlBase::XmlErrorHandler object should do for error handler
  // Make a specialized handler for content and install.
  // Then parse the file and catch any exceptions
  //  that propogate out
  //
  int errorCode = 0;
  try
  {
    //    SAXPrintHandlers handler(encodingName, unRepFlags);
    XmltorootHandler contentHandler(rootFile);
    xmlBase::XmlErrorHandler  errorHandler;
    parser->setDocumentHandler(&contentHandler);
    parser->setErrorHandler(&errorHandler);
    parser->parse(xmlFile);
    errorCount = parser->getErrorCount();
  }
  catch (const OutOfMemoryException&)
  {
    std::cerr << "OutOfMemoryException" 
                              << std::endl;
    errorCode = 5;
  }
  catch (const XMLException& toCatch)
  {
    char* charMsg = XMLString::transcode(toCatch.getMessage());
    std::cerr << "\nAn error occurred\n  Error: "
              << charMsg
              << "\n" << std::endl;
    errorCode = 4;
    XMLString::release(&charMsg);
  }
  if(errorCode) {
    XMLPlatformUtils::Terminate();
    return errorCode;
  }

  // Delete the parser itself.  Must be done prior to calling Terminate, below.
  delete parser;

  // And call the termination method
  XMLPlatformUtils::Terminate();

  if (errorCount > 0)
    return 4;
  else
    return 0;
}


