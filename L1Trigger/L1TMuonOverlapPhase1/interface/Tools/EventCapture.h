/*
 * EventCapture.h
 *
 *  Created on: Oct 23, 2019
 *      Author: kbunkow
 */

#ifndef OMTF_EVENTCAPTURE_H_
#define OMTF_EVENTCAPTURE_H_

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/CandidateSimMuonMatcher.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1I.h"
#include "TH2I.h"

class RPCGeometry;
class CSCGeometry;
class DTGeometry;

class EventCapture : public IOMTFEmulationObserver {
public:
  EventCapture(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
               CandidateSimMuonMatcher* candidateSimMuonMatcher);

  ~EventCapture() override;

  void beginRun(edm::EventSetup const& eventSetup) override;

  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>&,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const std::vector<l1t::RegionalMuonCand>& candMuons) override;

  void observeEventBegin(const edm::Event& event) override;

  void observeEventEnd(const edm::Event& event,
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override;

  void stubsSimHitsMatching(const edm::Event& iEvent);

  void endJob() override;

private:
  edm::InputTag simTrackInputTag;
  const OMTFConfiguration* omtfConfig;

  CandidateSimMuonMatcher* candidateSimMuonMatcher;

  std::vector<edm::Ptr<SimTrack> > simMuons;

  std::vector<std::shared_ptr<OMTFinput> > inputInProcs;
  std::vector<AlgoMuons> algoMuonsInProcs;
  std::vector<AlgoMuons> gbCandidatesInProcs;

  edm::InputTag rpcSimHitsInputTag;
  edm::InputTag cscSimHitsInputTag;
  edm::InputTag  dtSimHitsInputTag;

  edm::InputTag rpcDigiSimLinkInputTag;
  edm::InputTag cscStripDigiSimLinksInputTag;
  edm::InputTag dtDigiSimLinksInputTag;

  // pointers to the current geometry records
  unsigned long long _geom_cache_id = 0;
  edm::ESHandle<RPCGeometry> _georpc;
  edm::ESHandle<CSCGeometry> _geocsc;
  edm::ESHandle<DTGeometry> _geodt;

  TH2I* muonVsNotMuonStubs = nullptr;
  TH1I* muonStubsInLayers = nullptr;
  TH1I* notMuonStubsInLayers = nullptr;
};

#endif /* OMTF_EVENTCAPTURE_H_ */
