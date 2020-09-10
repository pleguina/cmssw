#ifndef OMTFProducer_H
#define OMTFProducer_H

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFReconstruction.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFSorter.h"

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

#include "SimDataFormats/Track/interface/SimTrackContainer.h"


class L1TMuonBayesOmtfTrackProducer : public edm::EDProducer {
 public:
  L1TMuonBayesOmtfTrackProducer(const edm::ParameterSet&);

  ~L1TMuonBayesOmtfTrackProducer() override;

  void beginJob() override;

  void endJob() override;

  void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  
  edm::EDGetTokenT<edm::SimTrackContainer> inputTokenSimHit; //TODO remove

  MuStubsInputTokens muStubsInputTokens;

  OMTFReconstruction m_Reconstruction;

};

#endif
