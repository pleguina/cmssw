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

  //.....................rrrrrrrrccccdddddd
  //.....................765432109876543210
  firedLayersToQuality[0b000000110000000011] = 1;
  firedLayersToQuality[0b000000100000000011] = 1;
  firedLayersToQuality[0b000000010000000011] = 1;
  firedLayersToQuality[0b000000110000000001] = 1;
  firedLayersToQuality[0b000001000000001100] = 1;
  firedLayersToQuality[0b000011000000001100] = 1;
  firedLayersToQuality[0b000010000000001100] = 1;
  firedLayersToQuality[0b000011000000000100] = 1;
  firedLayersToQuality[0b000000011000000001] = 1;
  firedLayersToQuality[0b001000010000000001] = 1;

  firedLayersToQuality[0b000100000000110000] = 1;
  firedLayersToQuality[0b001100000000010000] = 1;

  firedLayersToQuality[0b010000110000000001] = 8;
  firedLayersToQuality[0b000000111110000001] = 8;
  firedLayersToQuality[0b000000001000000011] = 8;
  firedLayersToQuality[0b000000111000000001] = 8;
  firedLayersToQuality[0b000000101000000001] = 8;
  firedLayersToQuality[0b010000011000000001] = 8;
  firedLayersToQuality[0b010000100000000001] = 8;
  firedLayersToQuality[0b000000110100000001] = 8;
  firedLayersToQuality[0b000000100100000001] = 8;
  firedLayersToQuality[0b001000100000000001] = 8;
  firedLayersToQuality[0b010000010000000001] = 8;
  firedLayersToQuality[0b001000110000000001] = 8;
  firedLayersToQuality[0b001000110000000000] = 8;
  firedLayersToQuality[0b000000010100000001] = 8;
  firedLayersToQuality[0b000010100000000001] = 8;
  firedLayersToQuality[0b000000100010000001] = 8;
  firedLayersToQuality[0b001010010000000101] = 8;
  firedLayersToQuality[0b100000000000000011] = 8;
  firedLayersToQuality[0b011011000000000000] = 8;
  firedLayersToQuality[0b000010110000000001] = 8;
  firedLayersToQuality[0b001001110000000001] = 8;
  firedLayersToQuality[0b000010100000000101] = 8;
  firedLayersToQuality[0b000011110000000001] = 8;
  firedLayersToQuality[0b000011110000001101] = 8;
  firedLayersToQuality[0b000011100000000101] = 8;
  firedLayersToQuality[0b000011110000000101] = 8;
  firedLayersToQuality[0b000100000001110000] = 8;
  firedLayersToQuality[0b000001110000001101] = 8;
  firedLayersToQuality[0b000000110110000001] = 8;
  firedLayersToQuality[0b000001110000000001] = 8;
  firedLayersToQuality[0b001000010001000001] = 8;
  firedLayersToQuality[0b000001100000000101] = 8;
  firedLayersToQuality[0b000001100000000001] = 8;
  firedLayersToQuality[0b000001110000000101] = 8;
  firedLayersToQuality[0b001001110001000001] = 8;
  firedLayersToQuality[0b000010110000000101] = 8;
  firedLayersToQuality[0b000000010001000001] = 8;
  firedLayersToQuality[0b000000100110000001] = 8;
  firedLayersToQuality[0b001001100000001100] = 8;
  firedLayersToQuality[0b000001010000000001] = 8;
  firedLayersToQuality[0b000010100000000011] = 8;
  firedLayersToQuality[0b000000100001000001] = 8;
  firedLayersToQuality[0b001000110001000001] = 8;
  firedLayersToQuality[0b000000010010000001] = 8;
  firedLayersToQuality[0b000001010000000101] = 8;
  firedLayersToQuality[0b100000100110000001] = 8;
  firedLayersToQuality[0b000010010000000101] = 8;
  firedLayersToQuality[0b000000110010000001] = 8;
  firedLayersToQuality[0b000000000000110100] = 8;
  firedLayersToQuality[0b000000010000000101] = 8;
  firedLayersToQuality[0b000000110001000001] = 8;
  firedLayersToQuality[0b000000010000001100] = 8;
  firedLayersToQuality[0b000010110000001101] = 8;
  firedLayersToQuality[0b000011010000001101] = 8;
  firedLayersToQuality[0b000000100000010001] = 8;
  firedLayersToQuality[0b000000110000000101] = 8;
  firedLayersToQuality[0b000001100000000111] = 8;
  firedLayersToQuality[0b000000100000000101] = 8;
  firedLayersToQuality[0b010000010010000001] = 8;
  firedLayersToQuality[0b000001100000001101] = 8;
  firedLayersToQuality[0b000011100000000111] = 8;
  firedLayersToQuality[0b000000010110000001] = 8;
  firedLayersToQuality[0b000011110000000111] = 8;
  firedLayersToQuality[0b000000011100000000] = 8;
  firedLayersToQuality[0b001000010000000011] = 8;
  firedLayersToQuality[0b000001110000000011] = 8;
  firedLayersToQuality[0b000100000000110000] = 8;
  firedLayersToQuality[0b000111100000110100] = 8;
  firedLayersToQuality[0b010000010010000000] = 8;
  firedLayersToQuality[0b100000010100000000] = 8;
  firedLayersToQuality[0b001000100000000011] = 8;
  firedLayersToQuality[0b000011100000001101] = 8;
  firedLayersToQuality[0b100000011100000000] = 8;
  firedLayersToQuality[0b110000011110000001] = 8;
  //firedLayersToQuality[0b000000000000110011] = 8;
  //firedLayersToQuality[0b000000100110000011] = 8;
  //firedLayersToQuality[0b110000000100000000] = 8;
  //firedLayersToQuality[0b001011110001001101] = 8;
  //firedLayersToQuality[0b010000100001000011] = 8;
  //firedLayersToQuality[0b000001100000001100] = 8;
  //firedLayersToQuality[0b000001110001000011] = 8;
  //firedLayersToQuality[0b011000000010000000] = 8;
  //firedLayersToQuality[0b001000110100000011] = 8;
  //firedLayersToQuality[0b010001000011000000] = 8;
  //firedLayersToQuality[0b100000000110000000] = 8;
  //firedLayersToQuality[0b000000000000111100] = 8;
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

  //TODO un-comment when convertToOuputScalesPhase2 is implemented
  /*
  omtfProc->setOutpuConversionFunction([&](l1t::tftype mtfType, const AlgoMuons& gbCandidates) {
    return this->convertToOuputScalesPhase2(mtfType, gbCandidates);
  }); */
}

void OmtfEmulation::getQualityFromFiredLayers(FinalMuon& finalMuon) {
  auto it = firedLayersToQuality.find(finalMuon.getAlgoMuon()->getFiredLayerBits());
  if (it != firedLayersToQuality.end()) {
    finalMuon.setQuality(it->second);
  } else {
    finalMuon.setQuality(12);  //default value
  }
};

FinalMuons OmtfEmulation::convertToOuputScalesPhase2(unsigned int iProcessor,
                                                     l1t::tftype mtfType,
                                                     const AlgoMuons& gbCandidates) {
  FinalMuons finalMuons;
  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>*>(omtfProc.get());
  if (omtfProcGoldenPat) {
    finalMuons = omtfProcGoldenPat->convertToOuputScalesPhase1(
        iProcessor, mtfType, gbCandidates);  //temporary solution, TODO remove

    if (ptAssignment) {
      for (auto& finalMuon : finalMuons) {
        //TODO convert the pts to the GMT output scales
        finalMuon.setPt(finalMuon.getAlgoMuon()->getPtNNConstr());
        finalMuon.setPtUnconstr(finalMuon.getAlgoMuon()->getPtNNUnconstr());
        finalMuon.setSign(finalMuon.getAlgoMuon()->getChargeNNConstr() < 0 ? 1 : 0);
        //finalMuon.setQuality(finalMuon.getAlgoMuon()->getQualityNN());

        getQualityFromFiredLayers(finalMuon);
      }
    }
    //TODO add conversion of eta anf phi from gbCandidates to the GMT output scales
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

    //TODO remove the below conversions when the conversions in the convertToOuputScalesPhase2 are implemented.
    ///N.B. the below conversions are from phase-1 uGMT scales to the phase-2 GMT scales.
    //What is needed in the convertToOuputScalesPhase2 is conversion from OTMF internal scales to  the phase-2 GMT scales
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
