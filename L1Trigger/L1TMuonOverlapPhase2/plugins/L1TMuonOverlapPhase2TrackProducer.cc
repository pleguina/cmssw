#include "L1Trigger/L1TMuonOverlapPhase2/plugins/L1TMuonOverlapPhase2TrackProducer.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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
           consumes<L1Phase2MuDTThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTThPhase2")),
           //pairs are consumed only if dtPairPolicy is set in the python config
           (edmParameterSet.exists("dtPairPolicy") && edmParameterSet.getParameter<bool>("dtPairPolicy"))
               ? consumes<L1Phase2MuDTExtPhiThetaPairContainer>(
                     edmParameterSet.getParameter<edm::InputTag>("srcDTPairsPhase2"))
               : edm::EDGetTokenT<L1Phase2MuDTExtPhiThetaPairContainer>()}),
      omtfParamsEsToken(esConsumes<L1TMuonOverlapParams, L1TMuonOverlapParamsRcd, edm::Transition::BeginRun>()),
      muonGeometryTokens({esConsumes<RPCGeometry, MuonGeometryRecord, edm::Transition::BeginRun>(),
                          esConsumes<CSCGeometry, MuonGeometryRecord, edm::Transition::BeginRun>(),
                          esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>()}),
      //needed for pattern generation and RootDataDumper
      magneticFieldEsToken(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagatorEsToken(esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(
          edm::ESInputTag("", "SteppingHelixPropagatorAlong"))),
      omtfEmulation(edmParameterSet, muStubsInputTokens, muStubsPhase2InputTokens) {
  produces<l1t::SAMuonCollection>("constrCands");
  produces<l1t::SAMuonCollection>("unconstrCands");

  //it is needed for pattern generation and RootDataDumper
  if (edmParameterSet.exists("simTracksTag"))
    mayConsume<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("simTracksTag"));
  if (edmParameterSet.exists("simVertexesTag"))
    mayConsume<edm::SimVertexContainer>(edmParameterSet.getParameter<edm::InputTag>("simVertexesTag"));
  if (edmParameterSet.exists("trackingParticleTag"))
    mayConsume<TrackingParticleCollection>(edmParameterSet.getParameter<edm::InputTag>("trackingParticleTag"));

  //truth and pileup collections read by DataROOTDumper2AllInput
  if (edmParameterSet.exists("genParticleTag"))
    mayConsume<reco::GenParticleCollection>(edmParameterSet.getParameter<edm::InputTag>("genParticleTag"));
  if (edmParameterSet.exists("pileupInfoInputTag"))
    mayConsume<std::vector<PileupSummaryInfo>>(edmParameterSet.getParameter<edm::InputTag>("pileupInfoInputTag"));

  //the dumper falls back to the g4SimHits muon hits when no explicit SimHit tag is given
  mayConsume<edm::PSimHitContainer>(edmParameterSet.exists("rpcSimHitsInputTag")
                                        ? edmParameterSet.getParameter<edm::InputTag>("rpcSimHitsInputTag")
                                        : edm::InputTag("g4SimHits", "MuonRPCHits"));
  mayConsume<edm::PSimHitContainer>(edmParameterSet.exists("cscSimHitsInputTag")
                                        ? edmParameterSet.getParameter<edm::InputTag>("cscSimHitsInputTag")
                                        : edm::InputTag("g4SimHits", "MuonCSCHits"));
  mayConsume<edm::PSimHitContainer>(edmParameterSet.exists("dtSimHitsInputTag")
                                        ? edmParameterSet.getParameter<edm::InputTag>("dtSimHitsInputTag")
                                        : edm::InputTag("g4SimHits", "MuonDTHits"));

  if (edmParameterSet.exists("rpcDigiSimLinkInputTag"))
    mayConsume<edm::DetSetVector<RPCDigiSimLink>>(edmParameterSet.getParameter<edm::InputTag>("rpcDigiSimLinkInputTag"));
  if (edmParameterSet.exists("cscStripDigiSimLinksInputTag"))
    mayConsume<edm::DetSetVector<StripDigiSimLink>>(
        edmParameterSet.getParameter<edm::InputTag>("cscStripDigiSimLinksInputTag"));
  if (edmParameterSet.exists("dtDigiSimLinksInputTag"))
    mayConsume<MuonDigiCollection<DTLayerId, DTDigiSimLink>>(
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

  auto outputCollections = omtfEmulation.run(iEvent, evSetup);

  iEvent.put(std::move(outputCollections.constrSaMuons), "constrCands");
  iEvent.put(std::move(outputCollections.unConstrSaMuons), "unconstrCands");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonOverlapPhase2TrackProducer);
