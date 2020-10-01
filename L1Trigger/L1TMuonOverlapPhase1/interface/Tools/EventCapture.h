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
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class EventCapture: public IOMTFEmulationObserver {
public:
  EventCapture(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig);

  virtual ~EventCapture();

  virtual void observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const std::shared_ptr<OMTFinput>&,
      const AlgoMuons& algoCandidates,
      const AlgoMuons& gbCandidates,
      const std::vector<l1t::RegionalMuonCand> & candMuons);

  virtual void observeEventBegin(const edm::Event& event);

  virtual void observeEventEnd(const edm::Event& event, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates);

  virtual void endJob();

private:
  edm::InputTag simTrackInputTag;
  const OMTFConfiguration* omtfConfig;

  std::vector<edm::Ptr< SimTrack > > simMuons;

  std::vector<std::shared_ptr<OMTFinput> > inputInProcs;
  std::vector<AlgoMuons> algoMuonsInProcs;
  std::vector<AlgoMuons> gbCandidatesInProcs;
};

#endif /* OMTF_EVENTCAPTURE_H_ */
