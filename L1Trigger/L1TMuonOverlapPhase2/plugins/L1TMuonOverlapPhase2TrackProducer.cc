#include "L1Trigger/L1TMuonOverlapPhase2/plugins/L1TMuonOverlapPhase2TrackProducer.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLink.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ProductRegistryHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1TMuonPhase2/interface/SAMuon.h"

#include <algorithm>
#include <iostream>
#include <memory>

L1TMuonOverlapPhase2TrackProducer::L1TMuonOverlapPhase2TrackProducer(const edm::ParameterSet& edmParameterSet)
    : muStubsInputTokens(
          {mayConsume<L1MuDTChambPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPh")),
           mayConsume<L1MuDTChambThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTTh")),
           consumes<CSCCorrelatedLCTDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcCSC")),
           consumes<RPCDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcRPC"))}),
      muStubsPhase2InputTokens(
          {consumes<L1Phase2MuDTPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPhPhase2")),
           consumes<L1Phase2MuDTThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTThPhase2"))}),
      omtfParamsEsToken(esConsumes<L1TMuonOverlapParams, L1TMuonOverlapParamsRcd, edm::Transition::BeginRun>()),
      muonGeometryTokens({esConsumes<RPCGeometry, MuonGeometryRecord, edm::Transition::BeginRun>(),
                          esConsumes<CSCGeometry, MuonGeometryRecord, edm::Transition::BeginRun>(),
                          esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>()}),
      //needed for pattern generation and RootDataDumper
      magneticFieldEsToken(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagatorEsToken(esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(
          edm::ESInputTag("", "SteppingHelixPropagatorAlong"))),
      omtfEmulation(edmParameterSet, muStubsInputTokens, muStubsPhase2InputTokens) {
  produces<l1t::RegionalMuonCandBxCollection>("OMTF");
  produces<l1t::SAMuonCollection>("OMTF");  //phase-2 collection

  //it is needed for pattern generation and RootDataDumper
  if (edmParameterSet.exists("simTracksTag"))
    mayConsume<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("simTracksTag"));
  if (edmParameterSet.exists("simVertexesTag"))
    mayConsume<edm::SimVertexContainer>(edmParameterSet.getParameter<edm::InputTag>("simVertexesTag"));
  if (edmParameterSet.exists("trackingParticleTag"))
    mayConsume<TrackingParticleCollection>(edmParameterSet.getParameter<edm::InputTag>("trackingParticleTag"));

  if (edmParameterSet.exists("genParticleTag"))
    mayConsume<reco::GenParticleCollection>(edmParameterSet.getParameter<edm::InputTag>("genParticleTag"));

  if (edmParameterSet.exists("rpcSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edmParameterSet.getParameter<edm::InputTag>("rpcSimHitsInputTag"));
  if (edmParameterSet.exists("cscSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edmParameterSet.getParameter<edm::InputTag>("cscSimHitsInputTag"));
  if (edmParameterSet.exists("dtSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edmParameterSet.getParameter<edm::InputTag>("dtSimHitsInputTag"));
  // Keep fallback registrations for configs that provide digi-sim links but not explicit SimHit tags.
  if (!edmParameterSet.exists("rpcSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "MuonRPCHits"));
  if (!edmParameterSet.exists("cscSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "MuonCSCHits"));
  if (!edmParameterSet.exists("dtSimHitsInputTag"))
    mayConsume<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "MuonDTHits"));

  if (edmParameterSet.exists("rpcDigiSimLinkInputTag"))
    mayConsume<edm::DetSetVector<RPCDigiSimLink> >(
        edmParameterSet.getParameter<edm::InputTag>("rpcDigiSimLinkInputTag"));
  if (edmParameterSet.exists("cscStripDigiSimLinksInputTag"))
    mayConsume<edm::DetSetVector<StripDigiSimLink> >(
        edmParameterSet.getParameter<edm::InputTag>("cscStripDigiSimLinksInputTag"));
  if (edmParameterSet.exists("dtDigiSimLinksInputTag"))
    mayConsume<MuonDigiCollection<DTLayerId, DTDigiSimLink> >(
        edmParameterSet.getParameter<edm::InputTag>("dtDigiSimLinksInputTag"));
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase2TrackProducer::beginJob() { omtfEmulation.beginJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase2TrackProducer::endJob() { omtfEmulation.endJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase2TrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
  omtfEmulation.beginRun(run, iSetup, omtfParamsEsToken, muonGeometryTokens, magneticFieldEsToken, propagatorEsToken);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapPhase2TrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup) {
  std::ostringstream str;

  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates = std::make_unique<l1t::RegionalMuonCandBxCollection>();

  std::unique_ptr<l1t::SAMuonCollection> saMuons = omtfEmulation.run(iEvent, evSetup, candidates);

  iEvent.put(std::move(saMuons), "OMTF");
  iEvent.put(std::move(candidates), "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonOverlapPhase2TrackProducer);
