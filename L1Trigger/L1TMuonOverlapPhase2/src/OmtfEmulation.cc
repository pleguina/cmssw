/*
 * OmtfEmulation.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include <memory>

#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfEmulation.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/PtAssignmentNNRegression.h"

#include "DataFormats/L1TMuonPhase2/interface/Constants.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>

OmtfEmulation::OmtfEmulation(const edm::ParameterSet& edmParameterSet,
                             MuStubsInputTokens& muStubsInputTokens,
                             MuStubsPhase2InputTokens& muStubsPhase2InputTokens)
    : OMTFReconstruction(edmParameterSet, muStubsInputTokens), muStubsPhase2InputTokens(muStubsPhase2InputTokens) {}

void OmtfEmulation::beginJob() {
  if (edmParameterSet.exists("usePhase2DTPrimitives") && edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
    inputMaker = std::make_unique<InputMakerPhase2>(edmParameterSet,
                                                    muStubsInputTokens,
                                                    muStubsPhase2InputTokens,
                                                    omtfConfig.get(),
                                                    std::make_unique<OmtfPhase2AngleConverter>());
  } else {
    inputMaker = std::make_unique<OMTFinputMaker>(
        edmParameterSet, muStubsInputTokens, omtfConfig.get(), std::make_unique<OmtfAngleConverter>());
  }
}

void OmtfEmulation::addObservers(const MuonGeometryTokens& muonGeometryTokens,
                                 const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>& magneticFieldEsToken,
                                 const edm::ESGetToken<Propagator, TrackingComponentsRecord>& propagatorEsToken) {
  if (observers.empty()) {  //assuring it is done only at the first run
    OMTFReconstruction::addObservers(muonGeometryTokens, magneticFieldEsToken, propagatorEsToken);
    /*    if(edmParameterSet.exists("patternsPtAssignment") && edmParameterSet.getParameter<bool>("patternsPtAssignment")) {
      //std::string rootFileName = edmParameterSet.getParameter<std::string>("dumpHitsFileName");
      .emplace_back(std::make_unique<PatternsPtAssignment>(edmParameterSet, omtfConfig.get(), omtfProcGoldenPat->getPatterns(), ""));
    }*/
  }

  //addObservers is called in OMTFReconstruction::beginRun after the omtfProc is constructed, therefore here we can used omtfProc
  if (edmParameterSet.exists("neuralNetworkFile") && !ptAssignment) {
    edm::LogImportant("OMTFReconstruction") << "constructing PtAssignmentNNRegression" << std::endl;
    std::string neuralNetworkFile = edmParameterSet.getParameter<edm::FileInPath>("neuralNetworkFile").fullPath();
    ptAssignment = std::make_unique<PtAssignmentNNRegression>(edmParameterSet, omtfConfig.get(), neuralNetworkFile);
  }

  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>*>(omtfProc.get());
  if (omtfProcGoldenPat) {
    omtfProcGoldenPat->setPtAssignment(ptAssignment.get());
    //omtfProcGoldenPat can be constructed from scratch each run, so ptAssignment is set herer every run
  }

  //otherwise the default convertToOuputScalesPhase1 is used
  if (ptAssignment) {
    omtfProc->setOutpuConversionFunction([&](l1t::tftype mtfType, const AlgoMuons& gbCandidates) {
      return this->convertToOuputScalesPhase2(mtfType, gbCandidates);
    });
  }
}

FinalMuons OmtfEmulation::convertToOuputScalesPhase2(l1t::tftype mtfType, const AlgoMuons& gbCandidates) {
  FinalMuons finalMuons;
  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>*>(omtfProc.get());
  if (omtfProcGoldenPat) {
    finalMuons = omtfProcGoldenPat->convertToOuputScalesPhase1(mtfType, gbCandidates);  //temporary solution, TODO remove
    if (ptAssignment) {
      for (auto& finalMuon : finalMuons) {
        finalMuon.setPt(finalMuon.getAlgoMuon()->getPtNNConstr());
        finalMuon.setPtUnconstr(finalMuon.getAlgoMuon()->getPtNNUnconstr());
        finalMuon.setSign(finalMuon.getAlgoMuon()->getChargeNNConstr() < 0 ? 1 : 0);
        finalMuon.setQuality(finalMuon.getAlgoMuon()->getQualityNN());
      }
    }
  }
  return finalMuons;
}

l1t::SAMuonCollection OmtfEmulation::getSAMuons(unsigned int iProcessor,
                                                l1t::tftype mtfType,
                                                FinalMuons& finalMuons,
                                                bool uncostrainedPt) {
  l1t::SAMuonCollection saMuons;

  for (auto& finalMuon : finalMuons) {
    unsigned int qual = finalMuon.getQuality();
    int charge = finalMuon.getSign();

    unsigned int pt = 0;
    if (!uncostrainedPt && finalMuon.getPt() > 0)
      pt = round(finalMuon.getPt() * 0.5 / Phase2L1GMT::LSBpt);  // Phase-1 LSB 0.5GeV
    if (uncostrainedPt && finalMuon.getPtUnconstr() > 0)
      pt = round(finalMuon.getPtUnconstr() * 1.0 / Phase2L1GMT::LSBpt);  // Phase-1 LSB 1.0GeV!!

    // BEWARE: THIS CONVERSION IS ONLY VALID FOR OMTF
    constexpr double p1phiLSB = 2 * M_PI / 576;
    // From the uGMTConfiguration of OMTF. OMTF send in local phi!!
    // all others correspond to 120 degree sectors = 192 in int-scale
    int globPhi = iProcessor * 192 + finalMuon.getPhi();
    // first processor starts at CMS phi = 15 degrees (24 in int)... Handle wrap-around with %. Add 576 to make sure the number is positive
    globPhi = (globPhi + 600) % 576;
    int phi = round(globPhi * p1phiLSB / Phase2L1GMT::LSBphi);             // Phase-1 LSB (2*pi/576)
    int eta = round(finalMuon.getEta() * 0.010875 / Phase2L1GMT::LSBeta);  // Phase-1 LSB 0.010875

    // FIXME: Below are not well defined in phase1 GMT
    // Using the version from Correlator for now
    int z0 = 0;  // No tracks info in Phase 1
    // Use 2 bits with LSB = 30cm for BMTF and 25cm for EMTF currently, but subjet to change
    int d0 = finalMuon.getHwD0();

    //Here do not use the word format to GT but use the word format expected by GMT
    /*
    int bstart = 0;
    wordtype word(0);
    bstart = wordconcat<wordtype>(word, bstart, 1, 1);
    bstart = wordconcat<wordtype>(word, bstart, charge, 1);
    bstart = wordconcat<wordtype>(word, bstart, pt, BITSPT);
    bstart = wordconcat<wordtype>(word, bstart, phi, BITSPHI);
    bstart = wordconcat<wordtype>(word, bstart, eta, BITSETA);
    //  bstart = wordconcat<wordtype>(word, bstart, z0, BITSSAZ0); NOT YET SUPPORTED BY GMT
    bstart = wordconcat<wordtype>(word, bstart, d0, BITSSAD0);
    bstart = wordconcat<wordtype>(
        word, bstart, qual, 8);  //FOR NOW 8 bits to be efficienct with Ghost busting. THIS IS ***NOT*** THE FINAL QUALITY
*/

    // Calculate Lorentz Vector
    math::PtEtaPhiMLorentzVector p4(pt * Phase2L1GMT::LSBpt, eta * Phase2L1GMT::LSBeta, phi * Phase2L1GMT::LSBphi, 0.0);
    l1t::SAMuon saMuon(p4, charge, pt, eta, phi, z0, d0, qual);
    saMuon.setTF(mtfType);
    //samuon.setWord(word);

    if (saMuon.hwPt() > 0) {
      saMuons.push_back(saMuon);
    }
  }

  return saMuons;
}

std::unique_ptr<l1t::SAMuonCollection> OmtfEmulation::run(
    const edm::Event& iEvent,
    const edm::EventSetup& evSetup,
    std::unique_ptr<l1t::RegionalMuonCandBxCollection>& candidates) {
  LogTrace("l1tOmtfEventPrint") << "\n" << __FUNCTION__ << ":" << __LINE__ << " iEvent " << iEvent.id().event() << endl;
  inputMaker->loadAndFilterDigis(iEvent);

  for (auto& obs : observers) {
    obs->observeEventBegin(iEvent);
  }

  std::unique_ptr<l1t::SAMuonCollection> saMuons = std::make_unique<l1t::SAMuonCollection>();
  candidates->setBXRange(bxMin, bxMax);

  ///The order is important: first put omtf_pos candidates, then omtf_neg.
  for (int bx = bxMin; bx <= bxMax; bx++) {
    for (unsigned int iProcessor = 0; iProcessor < omtfConfig->nProcessors(); ++iProcessor) {
      FinalMuons finalMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_pos, bx, inputMaker.get(), observers);

      l1t::SAMuonCollection procSAMuons = getSAMuons(iProcessor, l1t::tftype::omtf_pos, finalMuons, false);

      //fill outgoing collection
      for (auto& saMuon : procSAMuons) {
        saMuons->push_back(saMuon);
      }

      std::vector<l1t::RegionalMuonCand> candMuons =
          omtfProc->getRegionalMuonCands(iProcessor, l1t::tftype::omtf_pos, finalMuons);
      for (auto& candMuon : candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    for (unsigned int iProcessor = 0; iProcessor < omtfConfig->nProcessors(); ++iProcessor) {
      FinalMuons finalMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_neg, bx, inputMaker.get(), observers);

      l1t::SAMuonCollection procSAMuons = getSAMuons(iProcessor, l1t::tftype::omtf_neg, finalMuons, false);
      //fill outgoing collection
      for (auto& saMuon : procSAMuons) {
        saMuons->push_back(saMuon);
      }

      std::vector<l1t::RegionalMuonCand> candMuons =
          omtfProc->getRegionalMuonCands(iProcessor, l1t::tftype::omtf_neg, finalMuons);
      for (auto& candMuon : candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    //edm::LogInfo("OMTFReconstruction") <<"OMTF:  Number of candidates in BX="<<bx<<": "<<candidates->size(bx) << std::endl;;
  }

  LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ << endl;
  for (auto& obs : observers) {
    obs->observeEventEnd(iEvent, candidates);
  }

  return saMuons;
}
