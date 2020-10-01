#ifndef OMTFProducer_H
#define OMTFProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfEmulation.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

class L1TMuonOmtfPhase2TrackProducer : public edm::EDProducer {
 public:
  L1TMuonOmtfPhase2TrackProducer(const edm::ParameterSet&);

  ~L1TMuonOmtfPhase2TrackProducer() override;

  void beginJob() override;

  void endJob() override;

  void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  
  edm::EDGetTokenT<edm::SimTrackContainer> inputTokenSimHit; //TODO remove

  MuStubsInputTokens muStubsInputTokens;

  OmtfEmulation omtfEmulation;

};

#endif
