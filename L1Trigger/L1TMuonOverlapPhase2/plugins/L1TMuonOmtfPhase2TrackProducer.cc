#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ProductRegistryHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "L1Trigger/L1TMuonOverlapPhase2/plugins/L1TMuonOmtfPhase2TrackProducer.h"

#include <algorithm>
#include <iostream>
#include <memory>

L1TMuonOmtfPhase2TrackProducer::L1TMuonOmtfPhase2TrackProducer(const edm::ParameterSet& edmParameterSet) :
  muStubsInputTokens(
    {
      consumes<L1MuDTChambPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPh")),
      consumes<L1MuDTChambThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTTh")),
      consumes<CSCCorrelatedLCTDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcCSC")),
      consumes<RPCDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcRPC"))
    } ),
    omtfEmulation(edmParameterSet, muStubsInputTokens,
        consumes<L1Phase2MuDTPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPhPhase2")) )
{
  produces<l1t::RegionalMuonCandBxCollection >("OMTF");

  inputTokenSimHit = consumes<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("g4SimTrackSrc")); //TODO remove
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TMuonOmtfPhase2TrackProducer::~L1TMuonOmtfPhase2TrackProducer(){  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOmtfPhase2TrackProducer::beginJob(){

  omtfEmulation.beginJob();

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOmtfPhase2TrackProducer::endJob(){

  omtfEmulation.endJob();

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOmtfPhase2TrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup){

  omtfEmulation.beginRun(run, iSetup);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOmtfPhase2TrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup){

  std::ostringstream str;
  
  std::unique_ptr<l1t::RegionalMuonCandBxCollection > candidates = omtfEmulation.reconstruct(iEvent, evSetup);

  iEvent.put(std::move(candidates), "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonOmtfPhase2TrackProducer);
