#ifndef OMTFHitAnalyzer_H
#define OMTFHitAnalyzer_H

#include "xercesc/util/XercesDefs.hpp"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include <TH1.h>

class OMTFProcessor;
class OMTFConfiguration;
class OMTFConfigMaker;
class OMTFinputMaker;

class SimTrack;

class XMLConfigWriter;

namespace XERCES_CPP_NAMESPACE{
  class DOMElement;
  class DOMDocument;
  class DOMImplementation;
}

class OMTFHitAnalyzer : public edm::EDAnalyzer {
public:

	OMTFHitAnalyzer(const edm::ParameterSet & cfg);

  virtual ~OMTFHitAnalyzer();

  virtual void beginRun(edm::Run const& run, edm::EventSetup const& iSetup);

  virtual void beginJob();

  virtual void endJob();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);  

private:
  void configureProcesor(const OMTFConfiguration * omtfConfig,
      const L1TMuonOverlapParams* omtfPatterns, OMTFProcessor* omtfProc, unsigned int ptCode, int charge, unsigned int patNum);

  const SimTrack *findSimMuon(const edm::Event &ev, const edm::EventSetup &es, const SimTrack *previous=0);

  edm::ParameterSet theConfig;
  edm::InputTag g4SimTrackSrc;

  edm::EDGetTokenT<L1MuDTChambPhContainer> inputTokenDTPh;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDTTh;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCSC;
  edm::EDGetTokenT<RPCDigiCollection> inputTokenRPC;
  edm::EDGetTokenT<edm::SimTrackContainer> inputTokenSimHit;

  void writeMergedGPs();
  
  bool makeConnectionsMaps, makeGoldenPatterns, mergeXMLFiles;

  ///Original pdf width. read from configuration.
  unsigned int nPdfAddrBits;
  
  ///OMTF objects
  OMTFConfiguration *myOMTFConfig;
  OMTFinputMaker *myInputMaker;
  OMTFProcessor *myOMTF;
  ///
  xercesc::DOMElement *aTopElement;
  OMTFConfigMaker *myOMTFConfigMaker;
  XMLConfigWriter *myWriter;

  TH1I* ptDist;

}; 

#endif
