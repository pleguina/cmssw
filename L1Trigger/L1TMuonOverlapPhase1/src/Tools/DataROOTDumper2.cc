/*
 * DataROOTDumper2.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper2.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TFile.h"
#include "TTree.h"

DataROOTDumper2::DataROOTDumper2(const edm::ParameterSet& edmCfg,
                                 const OMTFConfiguration* omtfConfig,
                                 CandidateSimMuonMatcher* candidateSimMuonMatcher)
    : EmulationObserverBase(edmCfg, omtfConfig), candidateSimMuonMatcher(candidateSimMuonMatcher) {
  edm::LogVerbatim("l1tOmtfEventPrint") << " omtfConfig->nTestRefHits() " << omtfConfig->nTestRefHits()
                                        << " event.omtfGpResultsPdfSum.num_elements() " << endl;
  initializeTTree();

  if (edmCfg.exists("dumpKilledOmtfCands"))
    if (edmCfg.getParameter<bool>("dumpKilledOmtfCands"))
      dumpKilledOmtfCands = true;

  edm::LogVerbatim("l1tOmtfEventPrint") << " DataROOTDumper2 created. dumpKilledOmtfCands " << dumpKilledOmtfCands
                                        << std::endl;
}

DataROOTDumper2::~DataROOTDumper2() {}

void DataROOTDumper2::initializeTTree() {
  edm::Service<TFileService> fs;

  rootTree = fs->make<TTree>("OMTFHitsTree", "");

  rootTree->Branch("runNum", &omtfEvent.runNum);
  rootTree->Branch("lumiNum", &omtfEvent.lumiNum);
  rootTree->Branch("eventNum", &omtfEvent.eventNum);
  rootTree->Branch("eventNum64", &omtfEvent.eventNum64);
  rootTree->Branch("muonEvent", &omtfEvent.muonEvent);
  rootTree->Branch("muonBunchCrossing", &omtfEvent.muonBunchCrossing);
  rootTree->Branch("hasGenMatch", &omtfEvent.hasGenMatch);

  rootTree->Branch("muonPt", &omtfEvent.muonPt);
  rootTree->Branch("muonEta", &omtfEvent.muonEta);
  rootTree->Branch("muonPhi", &omtfEvent.muonPhi);
  rootTree->Branch("muonPropEta", &omtfEvent.muonPropEta);
  rootTree->Branch("muonPropPhi", &omtfEvent.muonPropPhi);
  rootTree->Branch("muonCharge", &omtfEvent.muonCharge);

  rootTree->Branch("muonDxy", &omtfEvent.muonDxy);
  rootTree->Branch("muonRho", &omtfEvent.muonRho);
  rootTree->Branch("parentPdgId", &omtfEvent.parentPdgId);
  rootTree->Branch("vertexEta", &omtfEvent.vertexEta);
  rootTree->Branch("vertexPhi", &omtfEvent.vertexPhi);

  rootTree->Branch("omtfPt", &omtfEvent.omtfPt);
  rootTree->Branch("omtfUPt", &omtfEvent.omtfUPt);
  rootTree->Branch("omtfEta", &omtfEvent.omtfEta);
  rootTree->Branch("omtfPhi", &omtfEvent.omtfPhi);
  rootTree->Branch("omtfCharge", &omtfEvent.omtfCharge);

  rootTree->Branch("omtfHwEta", &omtfEvent.omtfHwEta);

  rootTree->Branch("omtfProcessor", &omtfEvent.omtfProcessor);
  rootTree->Branch("omtfScore", &omtfEvent.omtfScore);
  rootTree->Branch("omtfQuality", &omtfEvent.omtfQuality);
  rootTree->Branch("omtfRefLayer", &omtfEvent.omtfRefLayer);
  rootTree->Branch("omtfRefHitNum", &omtfEvent.omtfRefHitNum);
  rootTree->Branch("omtfRefHitPhi", &omtfEvent.omtfRefHitPhi);
  // Explicit leaflist ("/B") required here: ROOT's TTree::Branch(name, T*) template overload
  // does not recognize plain "signed char" as a fundamental type (only "char"/Char_t triggers
  // its built-in scalar dispatch), so without a leaflist it falls through to object/class-based
  // branch creation and aborts with "Fatal Root Error: Illegal data type". "B" = 8-bit signed
  // integer (Char_t), matching the same pattern already used for reg_mtfType/reg_endcap in
  // DataROOTDumper2AllInput.cc.
  rootTree->Branch("omtfEndcap", &omtfEvent.omtfEndcap, "omtfEndcap/B");
  rootTree->Branch("omtfBx", &omtfEvent.omtfBx, "omtfBx/B");
  rootTree->Branch("candidateId", &omtfEvent.candidateId);

  rootTree->Branch("omtfFiredLayers", &omtfEvent.omtfFiredLayers);  //<<<<<<<<<<<<<<<<<<<<<<!!!!TODOO

  rootTree->Branch("killed", &omtfEvent.killed);

  rootTree->Branch("omtfHwPtConstr", &omtfEvent.omtfHwPtConstr);
  rootTree->Branch("omtfHwPtUnconstr", &omtfEvent.omtfHwPtUnconstr);
  rootTree->Branch("omtfPtNNConstr", &omtfEvent.omtfPtNNConstr);
  rootTree->Branch("omtfChargeNNConstr", &omtfEvent.omtfChargeNNConstr);
  rootTree->Branch("omtfDxyNn", &omtfEvent.omtfDxyNn);
  rootTree->Branch("omtfNnOutputs", &omtfEvent.omtfNnOutputs);
  rootTree->Branch("omtfPdfSumConstr", &omtfEvent.omtfPdfSumConstr);
  rootTree->Branch("omtfPdfSumUnconstr", &omtfEvent.omtfPdfSumUnconstr);
  rootTree->Branch("omtfFiredLayerCntConstr", &omtfEvent.omtfFiredLayerCntConstr);
  rootTree->Branch("omtfFiredLayerCntUnconstr", &omtfEvent.omtfFiredLayerCntUnconstr);
  rootTree->Branch("omtfFiredLayersConstr", &omtfEvent.omtfFiredLayersConstr);
  rootTree->Branch("omtfFiredLayersUnconstr", &omtfEvent.omtfFiredLayersUnconstr);
  rootTree->Branch("omtfPatternNumConstr", &omtfEvent.omtfPatternNumConstr);
  rootTree->Branch("omtfPatternNumUnconstr", &omtfEvent.omtfPatternNumUnconstr);
  rootTree->Branch("omtfDisplacedOnly", &omtfEvent.omtfDisplacedOnly);
  rootTree->Branch("omtfGmtPt", &omtfEvent.omtfGmtPt);
  rootTree->Branch("omtfGmtPtUnconstr", &omtfEvent.omtfGmtPtUnconstr);
  rootTree->Branch("omtfGmtEta", &omtfEvent.omtfGmtEta);
  rootTree->Branch("omtfGmtPhi", &omtfEvent.omtfGmtPhi);
  rootTree->Branch("layerC_pdfVal", &omtfEvent.layerC_pdfVal);
  rootTree->Branch("layerU_pdfVal", &omtfEvent.layerU_pdfVal);
  rootTree->Branch("layerC_deltaPhi", &omtfEvent.layerC_deltaPhi);
  rootTree->Branch("layerU_deltaPhi", &omtfEvent.layerU_deltaPhi);
  rootTree->Branch("layerC_pdfBin", &omtfEvent.layerC_pdfBin);
  rootTree->Branch("layerU_pdfBin", &omtfEvent.layerU_pdfBin);
  rootTree->Branch("layerC_valid", &omtfEvent.layerC_valid);
  rootTree->Branch("layerU_valid", &omtfEvent.layerU_valid);
  rootTree->Branch("layer_hasStub", &omtfEvent.layer_hasStub);

  rootTree->Branch("hits", &omtfEvent.hits);

  // Extended per-stub branches (same ordering and length as 'hits')
  rootTree->Branch("hits_phiBHw", &omtfEvent.hits_phiBHw);
  rootTree->Branch("hits_phiHw",  &omtfEvent.hits_phiHw);
  rootTree->Branch("hits_r",      &omtfEvent.hits_r);
  rootTree->Branch("hits_type",   &omtfEvent.hits_type);
  rootTree->Branch("hits_bx",     &omtfEvent.hits_bx);
  rootTree->Branch("hits_primitiveId", &omtfEvent.hits_primitiveId);
  rootTree->Branch("hits_layer",       &omtfEvent.hits_layer);
  rootTree->Branch("hits_quality",     &omtfEvent.hits_quality);
  rootTree->Branch("hits_etaHw",       &omtfEvent.hits_etaHw);
  rootTree->Branch("hits_valid",       &omtfEvent.hits_valid);
  rootTree->Branch("hits_detId",       &omtfEvent.hits_detId);

  rootTree->Branch("deltaEta", &omtfEvent.deltaEta);
  rootTree->Branch("deltaPhi", &omtfEvent.deltaPhi);

  ptGenPos = fs->make<TH1I>("ptGenPos", "ptGenPos, eta at vertex 0.8 - 1.24", 400, 0, 200);  //TODO
  ptGenNeg = fs->make<TH1I>("ptGenNeg", "ptGenNeg, eta at vertex 0.8 - 1.24", 400, 0, 200);
}

void DataROOTDumper2::observeProcesorEmulation(unsigned int iProcessor,
                                               l1t::tftype mtfType,
                                               const std::shared_ptr<OMTFinput>&,
                                               const AlgoMuons& algoCandidates,
                                               const AlgoMuons& gbCandidates,
                                               const FinalMuons& finalMuons) {}

void DataROOTDumper2::observeEventEnd(const edm::Event& iEvent, FinalMuons& finalMuons) {
  /*
  int muonCharge = 0;
  if (simMuon) {
    if (std::abs(simMuon->momentum().eta()) < 0.8 || std::abs(simMuon->momentum().eta()) > 1.24)
      return;

    muonCharge = (std::abs(simMuon->type()) == 13) ? simMuon->type() / -13 : 0;
    if (muonCharge > 0)
      ptGenPos->Fill(simMuon->momentum().pt());
    else
      ptGenNeg->Fill(simMuon->momentum().pt());
  }

  if (simMuon == nullptr || !omtfCand->isValid())  //no sim muon or empty candidate
    return;

  omtfEvent.muonPt = simMuon->momentum().pt();
  omtfEvent.muonEta = simMuon->momentum().eta();

  //TODO add cut on ete if needed
    if(std::abs(event.muonEta) < 0.8 || std::abs(event.muonEta) > 1.24)
    return;

  omtfEvent.muonPhi = simMuon->momentum().phi();
  omtfEvent.muonCharge = muonCharge;  //TODO
   */

  std::vector<MatchingResult> matchingResults = candidateSimMuonMatcher->getMatchingResults();
  LogTrace("l1tOmtfEventPrint") << "\nDataROOTDumper2::observeEventEnd matchingResults.size() "
                                << matchingResults.size() << std::endl;

  // Reset candidateId assignment for this event (see OmtfEvent::candidateId doc).
  nextCandidateIdForKey_.clear();

  //candidateSimMuonMatcher should use the  trackingParticles, because the simTracks are not stored for the pile-up events

  //for some events there are more than one matchingResults,
  //Usually at least one them has  genPt 0, which means no simMuon was matched, so candidate is ghost (or fake)
  //so better is to to drop such event, as it is not sure if the correct simMuon was matched to the candidate.
  //So we assume here that when the propagation is not used it is a single mu sample and this filter has sense
  //the propagation is used for multi-muon sample, so then this filter cannot be used
  //TODO add a flag to enable this filter? Disable it if not needed
  //in the single Mu sample, if there are two matchingResults, and the second has genPt 0, so it is easy to fitler it out when reading the dump.
  //so better would be to remove this condition, as it may be activated unintentionly
  for (auto& matchingResult : matchingResults) {
    omtfEvent = OmtfEvent();
    omtfEvent.runNum = iEvent.id().run();
    omtfEvent.lumiNum = iEvent.id().luminosityBlock();
    omtfEvent.eventNum = iEvent.id().event();
    omtfEvent.eventNum64 = static_cast<uint64_t>(iEvent.id().event());

    if (matchingResult.trackingParticle) {
      auto trackingParticle = matchingResult.trackingParticle;

      omtfEvent.hasGenMatch = true;

      if (matchingResult.result == MatchingResult::ResultType::propagationFailed)
        omtfEvent.muonEvent = -2;
      else
        omtfEvent.muonEvent = trackingParticle->eventId().event();
      omtfEvent.muonBunchCrossing = static_cast<int16_t>(trackingParticle->eventId().bunchCrossing());

      omtfEvent.muonPt = trackingParticle->pt();
      omtfEvent.muonEta = trackingParticle->momentum().eta();
      omtfEvent.muonPhi = trackingParticle->momentum().phi();
      omtfEvent.muonPropEta = matchingResult.propagatedEta;
      omtfEvent.muonPropPhi = matchingResult.propagatedPhi;
      //omtfEvent.muonCharge = (std::abs(trackingParticle->pdgId()) == 13) ? trackingParticle->pdgId() / -13 : 0;
      omtfEvent.muonCharge = trackingParticle->charge();

      if (trackingParticle->parentVertex().isNonnull()) {
        omtfEvent.muonDxy = trackingParticle->dxy();
        omtfEvent.muonRho = trackingParticle->parentVertex()->position().Rho();

        // Normally, parentVertex should have only one sourceTrack (in this context a pion or kaon)
        if (!trackingParticle->parentVertex()->sourceTracks().empty()) {
          omtfEvent.parentPdgId = trackingParticle->parentVertex()->sourceTracks()[0]->pdgId();
          LogTrace("l1MuonAnalyzerOmtf") << " DataROOTDumper2 parentTrackPdgId " << omtfEvent.parentPdgId << std::endl;
        }
      }

      omtfEvent.vertexPhi = matchingResult.vertexPhi;
      omtfEvent.vertexEta = matchingResult.vertexEta;

      omtfEvent.deltaEta = matchingResult.deltaEta;
      omtfEvent.deltaPhi = matchingResult.deltaPhi;

      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd trackingParticle: eventId "
                                    << trackingParticle->eventId().event() << " pdgId " << std::setw(3)
                                    << trackingParticle->pdgId() << " trackId "
                                    << trackingParticle->g4Tracks().at(0).trackId() << " pt " << std::setw(9)
                                    << trackingParticle->pt()  //<<" Beta "<<simMuon->momentum().Beta()
                                    << " eta " << std::setw(9) << trackingParticle->momentum().eta() << " phi "
                                    << std::setw(9) << trackingParticle->momentum().phi() << std::endl;

      if (std::abs(omtfEvent.muonEta) > 0.8 && std::abs(omtfEvent.muonEta) < 1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    } else if (matchingResult.simTrack) {
      auto simTrack = matchingResult.simTrack;

      omtfEvent.hasGenMatch = true;

      if (matchingResult.result == MatchingResult::ResultType::propagationFailed)
        omtfEvent.muonEvent = -2;
      else
        omtfEvent.muonEvent = simTrack->eventId().event();
      omtfEvent.muonBunchCrossing = static_cast<int16_t>(simTrack->eventId().bunchCrossing());

      omtfEvent.muonPt = simTrack->momentum().pt();
      omtfEvent.muonEta = simTrack->momentum().eta();
      omtfEvent.muonPhi = simTrack->momentum().phi();
      omtfEvent.muonPropEta = matchingResult.propagatedEta;
      omtfEvent.muonPropPhi = matchingResult.propagatedPhi;
      omtfEvent.muonCharge = simTrack->charge();

      /*if (!simTrack->noVertex() && matchingResult.simVertex) {
        const math::XYZTLorentzVectorD& vtxPos = matchingResult.simVertex->position();
        omtfEvent.muonDxy = (-vtxPos.X() * simTrack->momentum().py() + vtxPos.Y() * simTrack->momentum().px()) /
                            simTrack->momentum().pt();
        omtfEvent.muonRho = vtxPos.Rho();
      }*/

      omtfEvent.muonDxy = matchingResult.muonDxy;
      omtfEvent.muonRho = matchingResult.muonRho;

      omtfEvent.vertexPhi = matchingResult.vertexPhi;
      omtfEvent.vertexEta = matchingResult.vertexEta;

      omtfEvent.deltaEta = matchingResult.deltaEta;
      omtfEvent.deltaPhi = matchingResult.deltaPhi;

      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd simTrack: eventId "
                                    << simTrack->eventId().event() << " pdgId " << std::setw(3)
                                    << simTrack->type()  //<< " trackId " << simTrack->g4Tracks().at(0).trackId()
                                    << " pt " << std::setw(9)
                                    << simTrack->momentum().pt()  //<<" Beta "<<simMuon->momentum().Beta()
                                    << " eta " << std::setw(9) << simTrack->momentum().eta() << " phi " << std::setw(9)
                                    << simTrack->momentum().phi() << std::endl;

      if (std::abs(omtfEvent.muonEta) > 0.8 && std::abs(omtfEvent.muonEta) < 1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    } else {
      omtfEvent.hasGenMatch = false;
      omtfEvent.muonEvent = -1;
      omtfEvent.muonBunchCrossing = 0;

      omtfEvent.muonPt = 0;

      omtfEvent.muonEta = 0;
      omtfEvent.muonPhi = 0;

      omtfEvent.muonPropEta = 0;
      omtfEvent.muonPropPhi = 0;

      omtfEvent.muonCharge = 0;  //TODO

      omtfEvent.muonDxy = 0;
      omtfEvent.muonRho = 0;
    }

    auto addOmtfCand = [&](FinalMuonPtr muonCand) {
      //The constrained measurement is always defined for a valid candidate. In principle here only valid candidates should be, so this check is redundant.
      if (muonCand->getAlgoMuon()->getPdfSumConstr() > 0 && muonCand->getAlgoMuon()->getFiredLayerCntConstr() >= 3)
        omtfEvent.omtfPt = muonCand->getPtGev();
      else if (muonCand->getAlgoMuon()->getPtUnconstr() > 0)
        //if myCand->getPdfSumConstr() == 0, the myCand->getPtConstr() might not be 0, see the end of GhostBusterPreferRefDt::select
        //but hwPt=0 means empty candidate, hwPt=1 means pt=0,
        //but omtfPt = 0 means empty candidate
        //therefore here we set omtfPt=0.5 GeV, if the PtUnconstr > 0
        //N.B it is different than in the OMTFProcessor<GoldenPatternType>::convertToOuputScalesPhase1, where hwPt=1
        omtfEvent.omtfPt = 0.5;
      else
        omtfEvent.omtfPt = omtfConfig->hwPtToGev(0);

      //for candidate with no unconstrained measurement, hardware upt = 0
      //so then omtfEvent.omtfUPt is -1
      omtfEvent.omtfUPt = muonCand->getPtUnconstrGev();
      //omtfEvent.omtfEta = omtfConfig->hwEtaToEta(procMuon->getEtaHw());
      omtfEvent.omtfEta = muonCand->getEtaRad();
      omtfEvent.omtfPhi = muonCand->getPhiRad();
      omtfEvent.omtfCharge = muonCand->getAlgoMuon()->getChargeConstr();
      omtfEvent.omtfScore = muonCand->getAlgoMuon()->getPdfSum();

      omtfEvent.omtfHwEta = muonCand->getAlgoMuon()->getEtaHw();

      omtfEvent.omtfFiredLayers = muonCand->getAlgoMuon()->getFiredLayerBits();

      {
        const auto& algoMuon = *(muonCand->getAlgoMuon());
        const auto& gpC = algoMuon.getGpResultConstr();
        const auto& gpU = algoMuon.getGpResultUnconstr();
        omtfEvent.omtfHwPtConstr = algoMuon.getPtConstr();
        omtfEvent.omtfHwPtUnconstr = algoMuon.getPtUnconstr();
        omtfEvent.omtfPtNNConstr = algoMuon.getPtNNConstr();
        omtfEvent.omtfChargeNNConstr = algoMuon.getChargeNNConstr();
        omtfEvent.omtfDxyNn = algoMuon.getDxyNn();
        omtfEvent.omtfNnOutputs.assign(algoMuon.getNnOutputs().begin(), algoMuon.getNnOutputs().end());
        omtfEvent.omtfPdfSumConstr = gpC.getPdfSum();
        omtfEvent.omtfPdfSumUnconstr = gpU.getPdfSumUnconstr();
        omtfEvent.omtfFiredLayerCntConstr = gpC.getFiredLayerCnt();
        omtfEvent.omtfFiredLayerCntUnconstr = gpU.getFiredLayerCnt();
        omtfEvent.omtfFiredLayersConstr = gpC.getFiredLayerBits();
        omtfEvent.omtfFiredLayersUnconstr = gpU.getFiredLayerBits();
        omtfEvent.omtfPatternNumConstr = algoMuon.getPatternNumConstr();
        omtfEvent.omtfPatternNumUnconstr = algoMuon.getPatternNumUnconstr();
        //same condition as in OmtfProcessorPhase2::convertToGmtScalesPhase2
        omtfEvent.omtfDisplacedOnly = (algoMuon.getPdfSumConstr() == 0 || algoMuon.getFiredLayerCntConstr() < 3) &&
                                      algoMuon.getPtUnconstr() > 0;
        omtfEvent.omtfGmtPt = muonCand->getPtGmt();
        omtfEvent.omtfGmtPtUnconstr = muonCand->getPtUnconstrGmt();
        omtfEvent.omtfGmtEta = muonCand->getEtaGmt();
        omtfEvent.omtfGmtPhi = muonCand->getPhiGmt();

        const unsigned int nLayers = gpC.getStubResults().size();
        auto fillLayers = [nLayers](const GoldenPatternResult& gp,
                                    std::vector<float>& pdfVal,
                                    std::vector<short>& deltaPhi,
                                    std::vector<short>& pdfBin,
                                    std::vector<uint8_t>& valid) {
          pdfVal.assign(nLayers, 0);
          deltaPhi.assign(nLayers, 0);
          pdfBin.assign(nLayers, 0);
          valid.assign(nLayers, 0);
          for (unsigned int iLayer = 0; iLayer < nLayers && iLayer < gp.getStubResults().size(); ++iLayer) {
            const auto& stubResult = gp.getStubResults()[iLayer];
            pdfVal[iLayer] = stubResult.getPdfVal();
            deltaPhi[iLayer] = stubResult.getDeltaPhi();
            pdfBin[iLayer] = stubResult.getPdfBin();
            valid[iLayer] = stubResult.getValid();
          }
        };
        fillLayers(gpC, omtfEvent.layerC_pdfVal, omtfEvent.layerC_deltaPhi, omtfEvent.layerC_pdfBin, omtfEvent.layerC_valid);
        fillLayers(gpU, omtfEvent.layerU_pdfVal, omtfEvent.layerU_deltaPhi, omtfEvent.layerU_pdfBin, omtfEvent.layerU_valid);
        omtfEvent.layer_hasStub.assign(nLayers, 0);
        for (unsigned int iLayer = 0; iLayer < nLayers; ++iLayer)
          omtfEvent.layer_hasStub[iLayer] = gpC.getStubResults()[iLayer].getMuonStub() != nullptr;
      }
      omtfEvent.omtfRefLayer = muonCand->getAlgoMuon()->getRefLayer();
      omtfEvent.omtfRefHitNum = muonCand->getAlgoMuon()->getRefHitNumber();

      omtfEvent.hits.clear();
      omtfEvent.hits_phiBHw.clear();
      omtfEvent.hits_phiHw.clear();
      omtfEvent.hits_r.clear();
      omtfEvent.hits_type.clear();
      omtfEvent.hits_bx.clear();
      omtfEvent.hits_primitiveId.clear();
      omtfEvent.hits_layer.clear();
      omtfEvent.hits_quality.clear();
      omtfEvent.hits_etaHw.clear();
      omtfEvent.hits_valid.clear();
      omtfEvent.hits_detId.clear();

      //TODO choose, which gpResult should be dumped
      //auto& gpResult = procMuon->getGpResultConstr();
      auto& gpResult = (muonCand->getAlgoMuon()->getGpResultUnconstr().getPdfSumUnconstr() >
                        muonCand->getAlgoMuon()->getGpResultConstr().getPdfSum())
                           ? muonCand->getAlgoMuon()->getGpResultUnconstr()
                           : muonCand->getAlgoMuon()->getGpResultConstr();

      omtfEvent.omtfRefHitPhi = gpResult.getRefHitPhi();

      /*
        edm::LogVerbatim("l1tOmtfEventPrint")<<"DataROOTDumper2:;observeEventEnd muonPt "<<event.muonPt<<" muonCharge "<<event.muonCharge
            <<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer<<" omtfPtCont "<<event.omtfPtCont
            <<std::endl;  */

      for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
        auto& stubResult = gpResult.getStubResults()[iLogicLayer];

        //TODO it is to have the hit if it is below the quality cut
        /*if (omtfConfig->isBendingLayer(iLogicLayer) && !stubResult.getMuonStub()) {
          auto& stubResult = gpResult.getStubResults()[iLogicLayer-1];
        }*/

        if (stubResult.getMuonStub()) {  //&& stubResult.getValid() //TODO!!!!!!!!!!!!!!!!1
          OmtfEvent::Hit hit;
          hit.layer = iLogicLayer;
          hit.quality = stubResult.getMuonStub()->qualityHw;
          //hit.eta = stubResult.getMuonStub()->etaHw;  //replaced by deltaR
          hit.valid = stubResult.getValid();

          unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[muonCand->getAlgoMuon()->getRefLayer()];

          // Store raw phi displacement: hitPhi - refHitPhi (suitable for GNN edge features).
          // For bending layers (DT phiB pseudo-layers), phiDist carries phiBHw instead —
          // the bending angle is the only meaningful scalar on those layers.
          {
            int hitPhi = stubResult.getMuonStub()->phiHw;
            int phiRefHit = 0;
            if (gpResult.getStubResults()[refLayerLogicNum].getMuonStub())
              phiRefHit = gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->phiHw;
            hit.phiDist = hitPhi - phiRefHit;

            if (omtfConfig->isBendingLayer(iLogicLayer)) {
              hit.phiDist = stubResult.getMuonStub()->phiBHw;
            }
          }

          if (refLayerLogicNum == iLogicLayer)
            hit.deltaR = stubResult.getMuonStub()->r - 413;  //r of the ref hit - r of RB1in
          else
            hit.deltaR = stubResult.getMuonStub()->r -
                         (gpResult.getStubResults()[refLayerLogicNum].getMuonStub()
                              ? gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->r
                              : 413);

          LogTrace("l1tOmtfEventPrint")
              << " muonPt " << omtfEvent.muonPt << " omtfPt " << omtfEvent.omtfPt << " RefLayer "
              << int(omtfEvent.omtfRefLayer) << " layer " << int(hit.layer) << " PdfBin " << stubResult.getPdfBin()
              << " hit.phiDist " << hit.phiDist << " valid " << int(hit.valid) << " "  //<<" phiDist "<<phiDist
              << " hit.deltaR "
              << hit.deltaR
              << endl;

          if (hit.phiDist > 504 || hit.phiDist < -512) {
            LogTrace("l1tOmtfEventPrint")
                << " muonPt " << omtfEvent.muonPt << " omtfPt " << omtfEvent.omtfPt << " RefLayer "
                << (int)omtfEvent.omtfRefLayer << " layer " << int(hit.layer) << " hit.phiDist " << hit.phiDist
                << " valid " << stubResult.getValid() << " !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
          }

          //hit.etaHw is char, so we must limit the value being assigned
          //it char range is ok with valueP1Scale
          //for the phase2 scale something will have to be done TODO
          if (stubResult.getMuonStub()->etaHw > 127)
            hit.etaHw = 127;
          else if (stubResult.getMuonStub()->etaHw < -127)
            hit.etaHw = -127;
          else
            hit.etaHw = stubResult.getMuonStub()->etaHw;

          omtfEvent.hits.push_back(hit.rawData);

          // ---------- extended per-stub fields ----------
          omtfEvent.hits_phiBHw.push_back(static_cast<short>(stubResult.getMuonStub()->phiBHw));
          omtfEvent.hits_phiHw.push_back(static_cast<short>(stubResult.getMuonStub()->phiHw));
          omtfEvent.hits_r.push_back(static_cast<short>(stubResult.getMuonStub()->r));
          omtfEvent.hits_type.push_back(static_cast<signed char>(stubResult.getMuonStub()->type));
          {
            int bxVal = stubResult.getMuonStub()->bx;
            if (bxVal > 127) bxVal = 127;
            else if (bxVal < -127) bxVal = -127;
            omtfEvent.hits_bx.push_back(static_cast<signed char>(bxVal));
          }
          // Stable per-primitive identity, mirrors OMTFAllInputTree's reg_stub_* (item 2/3).
          omtfEvent.hits_primitiveId.push_back(stubResult.getMuonStub()->primitiveId);
          omtfEvent.hits_layer.push_back(static_cast<signed char>(iLogicLayer));
          omtfEvent.hits_quality.push_back(static_cast<signed char>(stubResult.getMuonStub()->qualityHw));
          omtfEvent.hits_etaHw.push_back(hit.etaHw);
          omtfEvent.hits_valid.push_back(stubResult.getValid() ? 1 : 0);
          omtfEvent.hits_detId.push_back(static_cast<uint32_t>(stubResult.getMuonStub()->detId));
          // ----------------------------------------------
        }
      }

      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd adding omtfCand : " << std::endl;

      LogTrace("l1tOmtfEventPrint") << " hwPt " << matchingResult.muonCand->getAlgoMuon()->getPtConstr() << " hwSign "
                                    << matchingResult.muonCand->getAlgoMuon()->getChargeConstr() << " hwQual "
                                    << matchingResult.muonCand->getQuality() << " hwEta " << std::setw(4)
                                    << matchingResult.muonCand->getAlgoMuon()->getEtaHw() << std::setw(4) << " hwPhi "
                                    << matchingResult.muonCand->getAlgoMuon()->getPhi() << "    eta " << std::setw(9)
                                    << matchingResult.muonCand->getEtaRad() << " isKilled "
                                    << matchingResult.muonCand->getAlgoMuon()->isKilled() << " tRefLayer "
                                    << matchingResult.muonCand->getAlgoMuon()->getRefLayer() << " RefHitNumber "
                                    << matchingResult.muonCand->getAlgoMuon()->getRefHitNumber() << std::endl;
    };

    if (matchingResult.muonCand) {
      //all candidates are dumped, whatever their quality and also the displaced-only ones (PtConstr == 0)
      omtfEvent.omtfQuality = matchingResult.muonCand->getQuality();
      omtfEvent.killed = false;
      omtfEvent.omtfProcessor = matchingResult.muonCand->getProcessor();

      omtfEvent.omtfEndcap = 0;
      if (matchingResult.muonCand->trackFinderType() == l1t::omtf_neg) {
        omtfEvent.omtfProcessor *= -1;
        omtfEvent.omtfEndcap = -1;
      } else if (matchingResult.muonCand->trackFinderType() == l1t::omtf_pos) {
        omtfEvent.omtfEndcap = 1;
      }
      // This offline dumper only ever produces the BX=0 emulation window (see OmtfEvent::omtfBx doc).
      omtfEvent.omtfBx = 0;
      omtfEvent.candidateId = nextCandidateIdForKey_[std::make_tuple(
          omtfEvent.omtfEndcap, omtfEvent.omtfProcessor, omtfEvent.omtfBx)]++;

      addOmtfCand(matchingResult.muonCand);
      rootTree->Fill();

      if (dumpKilledOmtfCands) {
        //the killed candidates have no FinalMuon, so it is built here as in OMTFProcessor::getFinalMuons
        auto& parentCand = matchingResult.muonCand;
        for (auto& killedAlgoMuon : parentCand->getAlgoMuon()->getKilledMuons()) {
          if (killedAlgoMuon->isKilled() == false) {
            edm::LogVerbatim("l1tOmtfEventPrint") << " killedCand->isKilled() == false !!!!!!!!";
          }
          auto killedCand = std::make_shared<FinalMuon>(killedAlgoMuon);
          killedCand->setPtGev(killedAlgoMuon->getPdfSumConstr() > 0 ? omtfConfig->hwPtToGev(killedAlgoMuon->getPtConstr())
                                                                     : 0.0);
          killedCand->setPtUnconstrGev(omtfConfig->hwPtToGev(killedAlgoMuon->getPtUnconstr()));
          float etaSign = (parentCand->trackFinderType() == l1t::omtf_neg) ? -1 : 1;
          killedCand->setEtaRad(etaSign * omtfConfig->hwEtaToEta(killedAlgoMuon->getEtaHw()));
          killedCand->setPhiRad(omtfConfig->procPhiOmtfToPhiRad(parentCand->getProcessor(), killedAlgoMuon->getPhi()));
          killedCand->setProcessor(parentCand->getProcessor());
          killedCand->setTrackFinderType(parentCand->trackFinderType());

          omtfEvent.omtfQuality = 0;
          omtfEvent.killed = true;
          addOmtfCand(killedCand);
          omtfEvent.candidateId = nextCandidateIdForKey_[std::make_tuple(
              omtfEvent.omtfEndcap, omtfEvent.omtfProcessor, omtfEvent.omtfBx)]++;
          rootTree->Fill();
        }
      }
    } else if (omtfEvent.muonPt > 0) {  //checking if there was a simMuon
      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd no matching omtfCand" << std::endl;

      omtfEvent.omtfPt = 0;
      omtfEvent.omtfUPt = 0;
      omtfEvent.omtfEta = 0;
      omtfEvent.omtfPhi = 0;
      omtfEvent.omtfCharge = 0;
      omtfEvent.omtfScore = 0;

      omtfEvent.omtfHwEta = 0;

      omtfEvent.omtfFiredLayers = 0;
      omtfEvent.omtfRefLayer = 0;
      omtfEvent.omtfRefHitNum = 0;
      omtfEvent.omtfProcessor = 10;
      omtfEvent.omtfEndcap = 0;
      omtfEvent.omtfBx = 0;
      omtfEvent.candidateId = nextCandidateIdForKey_[std::make_tuple(
          omtfEvent.omtfEndcap, omtfEvent.omtfProcessor, omtfEvent.omtfBx)]++;

      omtfEvent.omtfQuality = 0;
      omtfEvent.killed = false;

      omtfEvent.hits.clear();

      rootTree->Fill();
    }
  }
  evntCnt++;
}

void DataROOTDumper2::endJob() { edm::LogVerbatim("l1tOmtfEventPrint") << " evntCnt " << evntCnt << endl; }
