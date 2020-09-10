#include "L1Trigger/L1TMuonOverlapPhase1/plugins/L1TMuonBayesOmtfTrackProducer.h"

#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ProductRegistryHelper.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <algorithm>
#include <iostream>
#include <memory>

L1TMuonBayesOmtfTrackProducer::L1TMuonBayesOmtfTrackProducer(const edm::ParameterSet& edmParameterSet) :
  muStubsInputTokens(
    {
      consumes<L1MuDTChambPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPh")),
      consumes<L1MuDTChambThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTTh")),
      consumes<CSCCorrelatedLCTDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcCSC")),
      consumes<RPCDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcRPC"))
    } ),
  m_Reconstruction(edmParameterSet, muStubsInputTokens)
{
  produces<l1t::RegionalMuonCandBxCollection >("OMTF");

  inputTokenSimHit = consumes<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("g4SimTrackSrc")); //TODO remove
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TMuonBayesOmtfTrackProducer::~L1TMuonBayesOmtfTrackProducer(){
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtfTrackProducer::beginJob(){

  m_Reconstruction.beginJob();

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtfTrackProducer::endJob(){

  m_Reconstruction.endJob();

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtfTrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup){

  m_Reconstruction.beginRun(run, iSetup);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtfTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup){

  std::ostringstream str;
  
  std::unique_ptr<l1t::RegionalMuonCandBxCollection > candidates = m_Reconstruction.reconstruct(iEvent, evSetup);

  iEvent.put(std::move(candidates), "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonBayesOmtfTrackProducer);
