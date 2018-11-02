#ifndef L1TMuonOverlapTTMergerTrackProducer_H
#define L1TMuonOverlapTTMergerTrackProducer_H

#include "xercesc/util/XercesDefs.hpp"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFReconstruction.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

class L1TMuonOverlapParams;
class OMTFConfiguration;
class OMTFConfigMaker;
class XMLConfigWriter;



namespace XERCES_CPP_NAMESPACE{
  class DOMElement;
  class DOMDocument;
  class DOMImplementation;
}


class L1TMuonOverlapTTMergerTrackProducer : public edm::EDProducer {
 public:
  L1TMuonOverlapTTMergerTrackProducer(const edm::ParameterSet&);

  ~L1TMuonOverlapTTMergerTrackProducer() override;

  void beginJob() override;

  void endJob() override;

  void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

 private:

  edm::ParameterSet theConfig;
  
  edm::EDGetTokenT<L1MuDTChambPhContainer> inputTokenDTPh;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDTTh;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCSC;
  edm::EDGetTokenT<RPCDigiCollection> inputTokenRPC;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken;

  edm::EDGetTokenT<edm::SimTrackContainer> inputTokenSimHit; //TODO remove

  bool dumpResultToXML, dumpDetailedResultToXML;

  OMTFReconstruction m_Reconstruction;

};

#endif
