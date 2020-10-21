#include <L1Trigger/L1TMuonOverlapPhase1/plugins/L1TMuonOverlapPhase1TrackProducer.h>
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ProductRegistryHelper.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <algorithm>
#include <iostream>
#include <memory>

L1TMuonOverlapPhase1TrackProducer::L1TMuonOverlapPhase1TrackProducer(const edm::ParameterSet& edmParameterSet)
    : muStubsInputTokens(
          {consumes<L1MuDTChambPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPh")),
           consumes<L1MuDTChambThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTTh")),
           consumes<CSCCorrelatedLCTDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcCSC")),
           consumes<RPCDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcRPC"))}),
      m_Reconstruction(edmParameterSet, muStubsInputTokens) {
  produces<l1t::RegionalMuonCandBxCollection>("OMTF");

  inputTokenSimHit =
      consumes<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("g4SimTrackSrc"));  //TODO remove
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TMuonOverlapPhase1TrackProducer::~L1TMuonOverlapPhase1TrackProducer() {}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase1TrackProducer::beginJob() { m_Reconstruction.beginJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase1TrackProducer::endJob() { m_Reconstruction.endJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase1TrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
  m_Reconstruction.beginRun(run, iSetup);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase1TrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup) {
  std::ostringstream str;

  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates = m_Reconstruction.reconstruct(iEvent, evSetup);

  iEvent.put(std::move(candidates), "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonOverlapPhase1TrackProducer);
