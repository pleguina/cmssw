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

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Common/interface/Ptr.h"

#include <sstream>

EventCapture::EventCapture(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
    CandidateSimMuonMatcher* candidateSimMuonMatcher)
    : omtfConfig(omtfConfig),
      candidateSimMuonMatcher(candidateSimMuonMatcher),
      inputInProcs(omtfConfig->processorCnt()),
      algoMuonsInProcs(omtfConfig->processorCnt()),
      gbCandidatesInProcs(omtfConfig->processorCnt()) {
  //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" omtfConfig->nProcessors() "<<omtfConfig->nProcessors()<<std::endl;
  if (edmCfg.exists("g4SimTrackSrc"))
    simTrackInputTag = edmCfg.getParameter<edm::InputTag>("g4SimTrackSrc");
  else
    edm::LogImportant("OMTFReconstruction")
        << "EventCapture::EventCapture: no InputTag g4SimTrackSrc found" << std::endl;

  //rpcSimHitsInputTag = edmCfg.getParameter<edm::InputTag>("MuonRPCHits");
  //rpcSimHitsInputTag = edm::InputTag("g4SimHits", "MuonRPCHits");
  rpcSimHitsInputTag = edmCfg.getParameter<edm::InputTag>("rpcSimHitsInputTag");
  cscSimHitsInputTag = edmCfg.getParameter<edm::InputTag>("cscSimHitsInputTag");
   dtSimHitsInputTag = edmCfg.getParameter<edm::InputTag>("dtSimHitsInputTag");
}

EventCapture::~EventCapture() {
  // TODO Auto-generated destructor stub
}

void EventCapture::beginRun(edm::EventSetup const& eventSetup) {
  const MuonGeometryRecord& geom = eventSetup.get<MuonGeometryRecord>();
  unsigned long long geomid = geom.cacheIdentifier();
  if (_geom_cache_id != geomid) {
    geom.get(_georpc);
    geom.get(_geocsc);
    geom.get(_geodt);
    _geom_cache_id = geomid;
  }
}


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

  if(candidateSimMuonMatcher) {
    std::vector<MatchingResult>  matchingResults = candidateSimMuonMatcher->getMatchingResults();
    edm::LogVerbatim("l1tOmtfEventPrint")<<"matchingResults.size() "<<matchingResults.size()<<std::endl;

    for(auto& matchingResult : matchingResults) {
      if(matchingResult.muonCand && matchingResult.muonCand->hwQual() >= 12 && matchingResult.muonCand->hwPt() > 41 && matchingResult.genPt < 20) {
        dump = true;

        if(matchingResult.simTrack) {
          auto simMuon = matchingResult.simTrack;
          ostr<<"simMuon: eventId "<<simMuon->eventId().event()<<" pdgId "<<std::setw(3)<<simMuon->type()
                          <<" pt "<<std::setw(9)<<simMuon->momentum().pt() //<<" Beta "<<simMuon->momentum().Beta()
                          <<" eta "<<std::setw(9)<<simMuon->momentum().eta()<<" phi "<<std::setw(9)<<simMuon->momentum().phi()
                          <<std::endl;
        }
        else {
          ostr<<"no simMuon ";
        }
        ostr<<"matched to: "<<std::endl;
        auto finalCandidate = matchingResult.muonCand;
        ostr<< " hwPt " << finalCandidate->hwPt() << " hwSign " << finalCandidate->hwSign() << " hwQual "
                << finalCandidate->hwQual() << " hwEta " << std::setw(4) << finalCandidate->hwEta() << std::setw(4) << " hwPhi "
                << finalCandidate->hwPhi() << "    eta " << std::setw(9) << (finalCandidate->hwEta() * 0.010875) << " phi "
                << std::endl;

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

            auto globalPhiRad = omtfConfig->procHwPhiToGlobalPhi(stub->phiHw, OMTFinputMaker::getProcessorPhiZero(omtfConfig, iProc%6));
            ostrInput << (*stub) <<" globalPhiRad "<<globalPhiRad<<std::endl;
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

  stubsSimHitsMatching(iEvent);
}

void EventCapture::stubsSimHitsMatching(const edm::Event& iEvent) {
  edm::LogVerbatim("l1tOmtfEventPrint")<<"stubsSimHitsMatching ---------------" << std::endl;

  edm::Handle<edm::PSimHitContainer > rpcSimHitsHandle;
  iEvent.getByLabel(rpcSimHitsInputTag, rpcSimHitsHandle);
  edm::LogVerbatim("l1tOmtfEventPrint")<<"rpcSimHitsHandle: size: " << rpcSimHitsHandle->size()<<std::endl;

  edm::Handle<edm::PSimHitContainer > dtSimHitsHandle;
  iEvent.getByLabel(dtSimHitsInputTag, dtSimHitsHandle);
  edm::LogVerbatim("l1tOmtfEventPrint")<<std::endl<<"dtSimHitsHandle: size: " << dtSimHitsHandle->size()<<std::endl;

  edm::Handle<edm::PSimHitContainer > cscSimHitsHandle;
  iEvent.getByLabel(cscSimHitsInputTag, cscSimHitsHandle);
  edm::LogVerbatim("l1tOmtfEventPrint")<<std::endl<<"cscSimHitsHandle: size: " << cscSimHitsHandle->size()<<std::endl;

  for (unsigned int iProc = 0; iProc < gbCandidatesInProcs.size(); iProc++) {
    OmtfName board(iProc);

    edm::LogVerbatim("l1tOmtfEventPrint") << "gbCandidates " << std::endl;
    for (auto& gbCandidate : gbCandidatesInProcs[iProc])
      if (gbCandidate->isValid()) {
        edm::LogVerbatim("l1tOmtfEventPrint") << board.name() << " " << *gbCandidate << std::endl;
        edm::LogVerbatim("l1tOmtfEventPrint") << gbCandidate->getGpResult() << std::endl << std::endl;
        auto& gpResult = gbCandidate->getGpResult();
        for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
            auto& stub = gpResult.getStubResults()[iLogicLayer].getMuonStub();
            if (stub) {
              if(omtfConfig->isBendingLayer(iLogicLayer) )
                continue;

              if (gpResult.isLayerFired(iLogicLayer)) {

              }

              DetId stubDetId(stub->detId);
              if (stubDetId.det() != DetId::Muon) {
                edm::LogError("l1tOmtfEventPrint") << "!!!!!!!!!!!!!!!!!!!!!!!!  PROBLEM: hit in unknown Det, detID: " << stubDetId.det() << std::endl;
                continue;
              }

              auto stubGlobalPhi = omtfConfig->procHwPhiToGlobalPhi(stub->phiHw, OMTFinputMaker::getProcessorPhiZero(omtfConfig, iProc%6));
              edm::LogVerbatim("l1tOmtfEventPrint") << (*stub) <<"\nstubGlobalPhi "<<stubGlobalPhi<<std::endl;

              int matchedMuonHits = 0;
              int matchedNotMuonHits = 0;

              switch (stubDetId.subdetId()) {
                case MuonSubdetId::RPC: {
                  RPCDetId rpcDetId(stubDetId);


                  for(auto& simHit : *(rpcSimHitsHandle.product()) ) {
                    if(stubDetId.rawId() == simHit.detUnitId()) {
                      const RPCRoll* roll = _georpc->roll(rpcDetId);
                      auto strip = roll->strip(simHit.localPosition());
                      double simHitStripGlobalPhi = (roll->toGlobal(roll->centreOfStrip((int)strip))).phi();

                      if( abs(stubGlobalPhi - simHitStripGlobalPhi) < 0.02) {
                        if(abs(simHit.particleType()) == 13)
                          matchedMuonHits++;
                        else {
                          matchedNotMuonHits++;
                        }
                      }

                      edm::LogVerbatim("l1tOmtfEventPrint")
                                <<" simHitStripGlobalPhi "<<std::setw(10)<<simHitStripGlobalPhi
                                <<" strip "<<strip
                                <<" particleType: "<<simHit.particleType()
                                <<" event: "<<simHit.eventId().event()
                                <<" trackId "<<simHit.trackId()
                                <<" processType "<<simHit.processType()
                                <<" detUnitId "<<simHit.detUnitId()<<" "<<rpcDetId
                                //<<" phiAtEntry "<<simHit.phiAtEntry()
                                //<<" thetaAtEntry "<<simHit.thetaAtEntry()
                                //<<" localPosition: phi "<<simHit.localPosition().phi()<<" eta "<<simHit.localPosition().eta()
                                <<" localPosition: x "<<std::setw(10)<<simHit.localPosition().x()<<" y "<<std::setw(10)<<simHit.localPosition().y()
                                <<" timeOfFlight "<<simHit.timeOfFlight()
                                <<std::endl;
                    }
                  }
                  break;
                }//----------------------------------------------------------------------
                case MuonSubdetId::DT: {
                  //DTChamberId dt(stubDetId);
                  for(auto& simHit : *(dtSimHitsHandle.product()) ) {
                    const DTLayer* layer = _geodt->layer(DTLayerId(simHit.detUnitId()));
                    const DTChamber* chamber = layer->chamber();
                    if(stubDetId.rawId() == chamber->id().rawId() ) {
                      //auto strip = layer->geometry()->strip(simHit.localPosition());
                      auto simHitGlobalPoint = layer->toGlobal(simHit.localPosition());

                      edm::LogVerbatim("l1tOmtfEventPrint")
                                        <<" simHitGlobalPoint.phi "<<std::setw(10)<<simHitGlobalPoint.phi()
                                        //<<" strip "<<strip
                                        <<" particleType: "<<simHit.particleType()
                                        <<" event: "<<simHit.eventId().event()
                                        <<" trackId "<<simHit.trackId()
                                        <<" processType "<<simHit.processType()
                                        <<" detUnitId "<<simHit.detUnitId()<<" "<<layer->id()
                                        //<<" phiAtEntry "<<simHit.phiAtEntry()
                                        //<<" thetaAtEntry "<<simHit.thetaAtEntry()
                                        //<<" localPosition: phi "<<simHit.localPosition().phi()<<" eta "<<simHit.localPosition().eta()
                                        <<" localPosition: x "<<std::setw(10)<<simHit.localPosition().x()<<" y "<<std::setw(10)<<simHit.localPosition().y()
                                        <<" timeOfFlight "<<simHit.timeOfFlight()
                                        <<std::endl;
                    }
                  }
                  break;
                }//----------------------------------------------------------------------
                case MuonSubdetId::CSC: {
                  //CSCDetId csc(stubDetId);
                  for(auto& simHit : *(cscSimHitsHandle.product()) ) {
                    const CSCLayer* layer = _geocsc->layer(CSCDetId(simHit.detUnitId()));
                    auto chamber = layer->chamber();
                    if(stubDetId.rawId() == chamber->id().rawId() ) {
                      auto simHitStrip = layer->geometry()->strip(simHit.localPosition());
                      auto simHitGlobalPoint = layer->toGlobal(simHit.localPosition());
                      auto simHitStripGlobalPhi = layer->centerOfStrip(round(simHitStrip)).phi();

                      edm::LogVerbatim("l1tOmtfEventPrint")
                                <<" simHit: gloablPoint phi "<<simHitGlobalPoint.phi()
                                <<" stripGlobalPhi "<<simHitStripGlobalPhi.phi()
                                <<" strip "<<simHitStrip
                                <<" particleType: "<<simHit.particleType()
                                <<" event: "<<simHit.eventId().event()
                                <<" trackId "<<simHit.trackId()
                                <<" processType "<<simHit.processType()
                                <<" detUnitId "<<simHit.detUnitId()<<" "<<layer->id()
                                //<<" phiAtEntry "<<simHit.phiAtEntry()
                                //<<" thetaAtEntry "<<simHit.thetaAtEntry()
                                <<" timeOfFlight "<<simHit.timeOfFlight()
                                //<<" localPosition: phi "<<simHit.localPosition().phi()<<" eta "<<simHit.localPosition().eta()
                                <<" x "<<simHit.localPosition().x()<<" y "<<simHit.localPosition().y()

                                <<std::endl;
                    }
                  }

                  break;
                }
              }
              edm::LogVerbatim("l1tOmtfEventPrint")<<""<<std::endl;
            }
          }

      }
  }
}

void EventCapture::endJob() {}
