/*
 * EventCapture.cpp
 *
 *  Created on: Oct 23, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EventCapture.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfName.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include <sstream>

EventCapture::EventCapture(const edm::ParameterSet& edmCfg,
                           const OMTFConfiguration* omtfConfig,
                           CandidateSimMuonMatcher* candidateSimMuonMatcher)
    : omtfConfig(omtfConfig),
      candidateSimMuonMatcher(candidateSimMuonMatcher),
      inputInProcs(omtfConfig->processorCnt()),
      algoMuonsInProcs(omtfConfig->processorCnt()),
      gbCandidatesInProcs(omtfConfig->processorCnt()),
      stubsSimHitsMatcher(edmCfg, omtfConfig) {
  //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" omtfConfig->nProcessors() "<<omtfConfig->nProcessors()<<std::endl;
  if (edmCfg.exists("simTrackInputTag"))
    simTrackInputTag = edmCfg.getParameter<edm::InputTag>("simTrackInputTag");
  else
    edm::LogImportant("OMTFReconstruction")
        << "EventCapture::EventCapture: no InputTag simTrackInputTag found" << std::endl;
}

EventCapture::~EventCapture() {
  // TODO Auto-generated destructor stub
}

void EventCapture::beginRun(edm::EventSetup const& eventSetup) { stubsSimHitsMatcher.beginRun(eventSetup); }

void EventCapture::observeEventBegin(const edm::Event& event) {
  simMuons.clear();

  if (!simTrackInputTag.label().empty()) {
    edm::Handle<edm::SimTrackContainer> simTraksHandle;
    event.getByLabel(simTrackInputTag, simTraksHandle);

    for (unsigned int iSimTrack = 0; iSimTrack != simTraksHandle->size(); iSimTrack++) {
      if (abs((*simTraksHandle.product())[iSimTrack].type()) == 13)
        simMuons.emplace_back(simTraksHandle, iSimTrack);
    }
  }

  for (auto& input : inputInProcs)
    input.reset();

  for (auto& algoMuonsInProc : algoMuonsInProcs)
    algoMuonsInProc.clear();

  for (auto& gbCandidatesInProc : gbCandidatesInProcs)
    gbCandidatesInProc.clear();
}

void EventCapture::observeProcesorEmulation(unsigned int iProcessor,
                                            l1t::tftype mtfType,
                                            const std::shared_ptr<OMTFinput>& input,
                                            const AlgoMuons& algoCandidates,
                                            const AlgoMuons& gbCandidates,
                                            const std::vector<l1t::RegionalMuonCand>& candMuons) {
  unsigned int procIndx = omtfConfig->getProcIndx(iProcessor, mtfType);

  inputInProcs[procIndx] = input;

  algoMuonsInProcs[procIndx] = algoCandidates;
  gbCandidatesInProcs[procIndx] = gbCandidates;
}

void EventCapture::observeEventEnd(const edm::Event& iEvent,
                                   std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  std::ostringstream ostr;
  //filtering

  bool dump = false;

  if (candidateSimMuonMatcher) {
    std::vector<MatchingResult> matchingResults = candidateSimMuonMatcher->getMatchingResults();
    edm::LogVerbatim("l1tOmtfEventPrint") << "matchingResults.size() " << matchingResults.size() << std::endl;

    //candidateSimMuonMatcher should use the  trackingParticles, because the simTracks are not stored for the pile-up events
    for (auto& matchingResult : matchingResults) {
      if (matchingResult.muonCand && matchingResult.muonCand->hwQual() >= 12 &&
          matchingResult.muonCand->hwPt() > 38) {  //&& matchingResult.genPt < 20
        dump = true;

        bool runStubsSimHitsMatcher = false;
        if (matchingResult.trackingParticle) {
          auto trackingParticle = matchingResult.trackingParticle;
          ostr << "trackingParticle: eventId " << trackingParticle->eventId().event() << " pdgId " << std::setw(3)
               << trackingParticle->pdgId() << " trackId " << trackingParticle->g4Tracks().at(0).trackId() << " pt "
               << std::setw(9) << trackingParticle->pt()  //<<" Beta "<<simMuon->momentum().Beta()
               << " eta " << std::setw(9) << trackingParticle->momentum().eta() << " phi " << std::setw(9)
               << trackingParticle->momentum().phi() << std::endl;
        } else {
          ostr << "no simMuon ";
          runStubsSimHitsMatcher = true;
        }
        ostr << "matched to: " << std::endl;
        auto finalCandidate = matchingResult.muonCand;
        ostr << " hwPt " << finalCandidate->hwPt() << " hwSign " << finalCandidate->hwSign() << " hwQual "
             << finalCandidate->hwQual() << " hwEta " << std::setw(4) << finalCandidate->hwEta() << std::setw(4)
             << " hwPhi " << finalCandidate->hwPhi() << "    eta " << std::setw(9)
             << (finalCandidate->hwEta() * 0.010875) << " phi " << std::endl;

        if (runStubsSimHitsMatcher)
          stubsSimHitsMatcher.match(iEvent, matchingResult.muonCand, matchingResult.procMuon, ostr);
      }
    }
  }

  /*
  bool wasSimMuInOmtfPos = false;
  bool wasSimMuInOmtfNeg = false;
  for(auto& simMuon : simMuons) {
    if( simMuon->eventId().event() == 0 && abs(simMuon->momentum().eta() ) > 0.82 && abs(simMuon->momentum().eta() ) < 1.24 && simMuon->momentum().pt() > 5) {
      ostr<<"SimMuon: eventId "<<simMuon->eventId().event()<<" pdgId "<<std::setw(3)<<simMuon->type()
                  <<" pt "<<std::setw(9)<<simMuon->momentum().pt() //<<" Beta "<<simMuon->momentum().Beta()
                  <<" eta "<<std::setw(9)<<simMuon->momentum().eta()<<" phi "<<std::setw(9)<<simMuon->momentum().phi()
                  <<std::endl;

      if(simMuon->momentum().eta() > 0)
        wasSimMuInOmtfPos = true;
      else
        wasSimMuInOmtfNeg = true;
    }
  }


  bool wasCandInNeg = false;
  bool wasCandInPos = false;


  for(auto& finalCandidate : *finalCandidates) {
    if(finalCandidate.trackFinderType() == l1t::tftype::omtf_neg && finalCandidate.hwQual() >= 12 && finalCandidate.hwPt() > 41)
      wasCandInNeg = true;

    if(finalCandidate.trackFinderType() == l1t::tftype::omtf_pos && finalCandidate.hwQual() >= 12 && finalCandidate.hwPt() > 41)
      wasCandInPos = true;
  }


  if( (wasSimMuInOmtfNeg && !wasCandInNeg) )
    dump = true;

  if( (wasSimMuInOmtfPos && !wasCandInPos) )
    dump = true;

  */

  /*  dump = true; ///TODO if presetn then dumps all events!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(!dump)
    return;

  bool dump = false;
  for (auto& finalCandidate : *finalCandidates) {
    if (finalCandidate.hwPt() < 41) {  //  finalCandidate.hwQual() >= 1  41
      dump = true;
    }
  }*/

  if (!dump)
    return;

  ///printing

  edm::LogVerbatim("l1tOmtfEventPrint") << "##################### EventCapture::observeEventEnd - dump of event "
                                        << iEvent.id() << " #####################################################"
                                        << std::endl;

  edm::LogVerbatim("l1tOmtfEventPrint") << ostr.str() << endl;  //printing sim muons

  edm::LogVerbatim("l1tOmtfEventPrint") << "finalCandidates " << std::endl;
  for (auto& finalCandidate : *finalCandidates) {
    int globHwPhi = (finalCandidate.processor()) * 96 + finalCandidate.hwPhi();
    // first processor starts at CMS phi = 15 degrees (24 in int)... Handle wrap-around with %. Add 576 to make sure the number is positive
    globHwPhi = (globHwPhi + 600) % 576;

    double globalPhi = globHwPhi * 2. * M_PI / 576;
    if (globalPhi > M_PI)
      globalPhi = globalPhi - (2. * M_PI);

    int layerHits = (int)finalCandidate.trackAddress().at(0);
    std::bitset<18> layerHitBits(layerHits);

    edm::LogVerbatim("l1tOmtfEventPrint")
        << " hwPt " << finalCandidate.hwPt() << " hwSign " << finalCandidate.hwSign() << " hwQual "
        << finalCandidate.hwQual() << " hwEta " << std::setw(4) << finalCandidate.hwEta() << std::setw(4) << " hwPhi "
        << finalCandidate.hwPhi() << "    eta " << std::setw(9) << (finalCandidate.hwEta() * 0.010875) << " phi "
        << std::setw(9) << globalPhi << " " << layerHitBits << " processor "
        << OmtfName(finalCandidate.processor(), finalCandidate.trackFinderType()) << std::endl;

    for (auto& trackAddr : finalCandidate.trackAddress()) {
      if (trackAddr.first >= 10)
        edm::LogVerbatim("l1tOmtfEventPrint") << "trackAddr first " << trackAddr.first << " second " << trackAddr.second
                                              << " ptGeV " << omtfConfig->hwPtToGev(trackAddr.second);
    }
  }
  edm::LogVerbatim("l1tOmtfEventPrint") << std::endl;

  for (unsigned int iProc = 0; iProc < inputInProcs.size(); iProc++) {
    OmtfName board(iProc);

    std::ostringstream ostrInput;
    if (inputInProcs[iProc]) {
      auto& omtfInput = *inputInProcs[iProc];
      int layersWithStubs = 0;
      for (auto& layer : omtfInput.getMuonStubs()) {
        for (auto& stub : layer) {
          bool layerFired = false;
          if (stub && (stub->type != MuonStub::Type::EMPTY)) {
            layerFired = true;

            auto globalPhiRad = omtfConfig->procHwPhiToGlobalPhi(
                stub->phiHw, OMTFinputMaker::getProcessorPhiZero(omtfConfig, iProc % 6));
            ostrInput << (*stub) << " globalPhiRad " << globalPhiRad << std::endl;
          }
          if (layerFired)
            layersWithStubs++;
        }
      }

      if (layersWithStubs != 0) {
        edm::LogVerbatim("l1tOmtfEventPrint") << "\niProcessor " << iProc << " " << board.name()
                                              << " **************************************************" << std::endl;
        edm::LogVerbatim("l1tOmtfEventPrint") << ostrInput.str() << std::endl;
      }

      if (layersWithStubs < 2)
        continue;

      edm::LogVerbatim("l1tOmtfEventPrint") << *inputInProcs[iProc] << std::endl;

      edm::LogVerbatim("l1tOmtfEventPrint") << "algoMuons " << std::endl;
      for (auto& algoMuon : algoMuonsInProcs[iProc]) {
        if (algoMuon->isValid()) {
          edm::LogVerbatim("l1tOmtfEventPrint") << board.name() << " " << *algoMuon << std::endl;
          edm::LogVerbatim("l1tOmtfEventPrint") << algoMuon->getGpResult() << std::endl << std::endl;
        }
      }

      edm::LogVerbatim("l1tOmtfEventPrint") << "gbCandidates " << std::endl;
      for (auto& gbCandidate : gbCandidatesInProcs[iProc])
        if (gbCandidate->isValid())
          edm::LogVerbatim("l1tOmtfEventPrint") << board.name() << " " << *gbCandidate << std::endl;

      edm::LogVerbatim("l1tOmtfEventPrint") << std::endl;
    }
  }

  edm::LogVerbatim("l1tOmtfEventPrint") << std::endl;
}

void EventCapture::endJob() { stubsSimHitsMatcher.endJob(); }
