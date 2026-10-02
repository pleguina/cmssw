/*
 * InputMakerPhase2.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfName.h"

#include <iostream>
#include <map>
#include <tuple>

/////////////////////////////////////
void DtPhase2DigiToStubsConverter::loadDigis(const edm::Event& event) {
  event.getByToken(inputTokenDtPh, dtPhDigis);
  event.getByToken(inputTokenDtTh, dtThDigis);
  if (dtPairPolicy)
    event.getByToken(inputTokenDtPairs, dtPairs);
  bunchCrossing = event.bunchCrossing();
}

void DtPhase2DigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                             unsigned int iProcessor,
                                             l1t::tftype procTyp,
                                             int bxFrom,
                                             int bxTo,
                                             std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers,
                                             XmlIOCache& xmlCache) {
  boost::property_tree::ptree procDataTree;

  std::map<std::string, boost::property_tree::ptree> chamberTrees;

  if (dtPairPolicy) {
    //slot = position of the pair among the pairs of the same chamber and bx, in the order of the pair collection
    //(the order set by the producer, i.e. the order in which the BMT would send the stubs); empty slots are counted
    std::map<std::tuple<int, int, int, int>, int> slotCounters;

    for (const auto& pair : dtPairs->getContainer()) {
      const auto& phi = pair.phiDigi();
      const auto& th = pair.thetaDigi();
      constexpr int emptyDigiBx = -100;
      bool hasPhi = phi.bxNum() != emptyDigiBx;
      bool hasTheta = th.bxNum() != emptyDigiBx;
      if (!hasPhi && !hasTheta)
        continue;

      int wh = hasPhi ? phi.whNum() : th.whNum();
      int sc = hasPhi ? phi.scNum() : th.scNum();
      int st = hasPhi ? phi.stNum() : th.stNum();
      int bx = hasPhi ? phi.bxNum() : th.bxNum();

      int slot = slotCounters[std::make_tuple(wh, sc, st, bx)]++;

      DTChamberId detid(wh, st, sc + 1);
      if (!acceptDigi(detid, iProcessor, procTyp))
        continue;
      if (bx - config->dtBxShift() < bxFrom || bx - config->dtBxShift() > bxTo)
        continue;

      MuonStub addedStub;
      bool stubAdded = false;
      if (hasPhi)  //a theta without phi is only dumped, the OMTF algorithm uses the phi stubs
        stubAdded = addDTpair(muonStubsInLayers, pair, iProcessor, procTyp, &addedStub);

      std::ostringstream chamberName;
      chamberName << "DT_Wh" << wh << "_St_" << st << "_Se_" << sc + 1;
      auto& dtChamberTree = chamberTrees[chamberName.str()];
      auto& dtP2PairDigi = dtChamberTree.add_child("dtP2PairDigi", boost::property_tree::ptree());
      dtP2PairDigi.add("<xmlattr>.whNum", wh);
      dtP2PairDigi.add("<xmlattr>.scNum", sc);
      dtP2PairDigi.add("<xmlattr>.stNum", st);
      dtP2PairDigi.add("<xmlattr>.bx", bx);
      dtP2PairDigi.add("<xmlattr>.slot", slot);
      dtP2PairDigi.add("<xmlattr>.pairQuality", pair.quality());
      dtP2PairDigi.add("<xmlattr>.hasPhi", hasPhi);
      dtP2PairDigi.add("<xmlattr>.hasTheta", hasTheta);
      //the stub the emulator made from this pair (the expected output of the firmware DT pair interface)
      dtP2PairDigi.add("<xmlattr>.stubValid", stubAdded);
      if (stubAdded) {
        dtP2PairDigi.add("<xmlattr>.stubPhi", addedStub.phiHw);
        dtP2PairDigi.add("<xmlattr>.stubPhiB", addedStub.phiBHw);
        dtP2PairDigi.add("<xmlattr>.stubQuality", addedStub.qualityHw);
        dtP2PairDigi.add("<xmlattr>.stubEta", addedStub.etaHw);
      }
      if (hasPhi) {
        dtP2PairDigi.add("<xmlattr>.slNum", phi.slNum());
        dtP2PairDigi.add("<xmlattr>.quality", phi.quality());
        dtP2PairDigi.add("<xmlattr>.rpcFlag", phi.rpcFlag());
        dtP2PairDigi.add("<xmlattr>.phi", phi.phi());
        dtP2PairDigi.add("<xmlattr>.phiBend", phi.phiBend());
        dtP2PairDigi.add("<xmlattr>.t0", phi.t0());
      }
      if (hasTheta) {
        dtP2PairDigi.add("<xmlattr>.thQuality", th.quality());
        dtP2PairDigi.add("<xmlattr>.k", th.k());
        dtP2PairDigi.add("<xmlattr>.z", th.z());
        dtP2PairDigi.add("<xmlattr>.thT0", th.t0());
      }
    }

    for (auto& chamberTree : chamberTrees) {
      chamberTree.second.add("<xmlattr>.name", chamberTree.first);
      procDataTree.add_child("dtChamber", chamberTree.second);
    }
    for (auto& obs : observers)
      obs->addProcesorData("dtData", procDataTree);
    return;
  }

  for (const auto& digiIt : *dtPhDigis->getContainer()) {
    DTChamberId detid(digiIt.whNum(), digiIt.stNum(), digiIt.scNum() + 1);

    LogTrace("l1tOmtfEventPrint") << " L1Phase2MuDTPhDigi:35: detid " << detid << " digi "
                                  << " whNum " << digiIt.whNum() << " scNum " << digiIt.scNum() << " stNum "
                                  << digiIt.stNum() << " slNum " << digiIt.slNum() << " quality " << digiIt.quality()
                                  << " rpcFlag " << digiIt.rpcFlag() << " phi " << digiIt.phi() << " phiBend "
                                  << digiIt.phiBend() << " digBx " << digiIt.bxNum() << " correctedBx "
                                  << digiIt.bxNum() - bunchCrossing << std::endl;

    ///Check it the data fits into given processor input range
    if (!acceptDigi(detid, iProcessor, procTyp))
      continue;

    // HACK for Phase-2  (DT TPs are centered in bX=20 in MC, and 13(?) in data)
    if (digiIt.bxNum() - config->dtBxShift() >= bxFrom && digiIt.bxNum() - config->dtBxShift() <= bxTo) {
      addDTphiDigi(muonStubsInLayers, digiIt, dtThDigis.product(), iProcessor, procTyp);

      std::ostringstream chamberName;
      //chamberName<<detid;
      chamberName << "DT_Wh" << digiIt.whNum() << "_St_" << digiIt.stNum() << "_Se_" << digiIt.scNum() + 1;
      auto& dtChamberTree = chamberTrees[chamberName.str()];
      auto& dtP2PhiDigi = dtChamberTree.add_child("dtP2PhiDigi", boost::property_tree::ptree());
      dtP2PhiDigi.add("<xmlattr>.whNum", digiIt.whNum());
      dtP2PhiDigi.add("<xmlattr>.scNum", digiIt.scNum());
      dtP2PhiDigi.add("<xmlattr>.stNum", digiIt.stNum());
      dtP2PhiDigi.add("<xmlattr>.slNum", digiIt.slNum());
      dtP2PhiDigi.add("<xmlattr>.quality", digiIt.quality());
      dtP2PhiDigi.add("<xmlattr>.rpcFlag", digiIt.rpcFlag());
      dtP2PhiDigi.add("<xmlattr>.phi", digiIt.phi());
      dtP2PhiDigi.add("<xmlattr>.phiBend", digiIt.phiBend());
    }
  }

  for (auto& thetaDigi : (*(dtThDigis->getContainer()))) {
    if (thetaDigi.bxNum() - config->dtBxShift() >= bxFrom && thetaDigi.bxNum() - config->dtBxShift() <= bxTo) {
      if (!mergePhiAndTheta) {
        addDTetaStubs(muonStubsInLayers, thetaDigi, iProcessor, procTyp);
      }

      std::ostringstream chamberName;
      //chamberName<<detid;
      chamberName << "DT_Wh" << thetaDigi.whNum() << "_St_" << thetaDigi.stNum() << "_Se_" << thetaDigi.scNum() + 1;
      auto& dtChamberTree = chamberTrees[chamberName.str()];
      auto& dtP2ThDigi = dtChamberTree.add_child("dtP2ThDigi", boost::property_tree::ptree());
      dtP2ThDigi.add("<xmlattr>.whNum", thetaDigi.whNum());
      dtP2ThDigi.add("<xmlattr>.scNum", thetaDigi.scNum());
      dtP2ThDigi.add("<xmlattr>.stNum", thetaDigi.stNum());
      dtP2ThDigi.add("<xmlattr>.quality", thetaDigi.quality());
      dtP2ThDigi.add("<xmlattr>.rpcFlag", thetaDigi.rpcFlag());
      dtP2ThDigi.add("<xmlattr>.k", thetaDigi.k());
      dtP2ThDigi.add("<xmlattr>.z", thetaDigi.z());
    }
  }

  for (auto& chamberTree : chamberTrees) {
    chamberTree.second.add("<xmlattr>.name", chamberTree.first);
    procDataTree.add_child("dtChamber", chamberTree.second);
  }

  for (auto& obs : observers)
    obs->addProcesorData("dtData", procDataTree);
}

//dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
void DtPhase2DigiToStubsConverterOmtf::addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                                                    const L1Phase2MuDTPhDigi& digi,
                                                    const L1Phase2MuDTThContainer* dtThDigis,
                                                    unsigned int iProcessor,
                                                    l1t::tftype procTyp) {
  DTChamberId detid(digi.whNum(), digi.stNum(), digi.scNum() + 1);
  //no config->dtBxShift() here, and also no shift inside getGlobalEtaPhase2
  int etaHw = angleConverter->getGlobalEtaPhase2(detid, dtThDigis, digi.bxNum());
  addPhiStub(muonStubsInLayers, digi, etaHw, iProcessor, procTyp);
}

bool DtPhase2DigiToStubsConverterOmtf::addDTpair(MuonStubPtrs2D& muonStubsInLayers,
                                                 const L1Phase2MuDTExtPhiThetaPair& pair,
                                                 unsigned int iProcessor,
                                                 l1t::tftype procTyp,
                                                 MuonStub* outStub) {
  const auto& digi = pair.phiDigi();
  DTChamberId detid(digi.whNum(), digi.stNum(), digi.scNum() + 1);
  int etaHw = angleConverter->getEtaFromThetaDigi(detid, pair.thetaDigi(), dtPairEtaWindowCheck, dtPairFixedEta);
  return addPhiStub(muonStubsInLayers, digi, etaHw, iProcessor, procTyp, dtPairWordPrecision, outStub);
}

bool DtPhase2DigiToStubsConverterOmtf::addPhiStub(MuonStubPtrs2D& muonStubsInLayers,
                                                  const L1Phase2MuDTPhDigi& digi,
                                                  int etaHw,
                                                  unsigned int iProcessor,
                                                  l1t::tftype procTyp,
                                                  bool wordPrecision,
                                                  MuonStub* outStub) {
  DTChamberId detid(digi.whNum(), digi.stNum(), digi.scNum() + 1);

  MuonStub stub;

  //converting the quality to the same encoding as in phase-1, as it is important for extrapolation
  if (digi.quality() >= 6)  // correlated stub
    stub.qualityHw = digi.quality() - 2;
  else if (digi.quality() >= 3) {  // 4 hit uncorrelated stub
    if (digi.slNum() == 3)
      stub.qualityHw = 3;
    else if (digi.slNum() == 1)
      stub.qualityHw = 2;
  } else {  //quality 1 (3 hits) or 2 (3+2 hits)
    if (digi.slNum() == 3)
      stub.qualityHw = 1;
    else if (digi.slNum() == 1)
      stub.qualityHw = 0;
  }

  //BMT->OMTF stub word: the word has no superlayer, so the uncorrelated 4-hit qualities map to one value
  if (wordPrecision && digi.quality() >= 3 && digi.quality() < 6)
    stub.qualityHw = 3;

  if (stub.qualityHw < config->getMinDtPhiQuality())
    return false;

  unsigned int hwNumber = config->getLayerNumber(detid.rawId());
  auto iter = config->getHwToLogicLayer().find(hwNumber);
  if (iter == config->getHwToLogicLayer().end())
    return false;

  unsigned int iLayer = iter->second;
  unsigned int iInput = OMTFinputMaker::getInputNumber(config, detid.rawId(), iProcessor, procTyp);
  //MuonStub& stub = muonStubsInLayers[iLayer][iInput];

  stub.type = MuonStub::DT_PHI_ETA;

  //the stub word carries phi with 13 bits, i.e. the 17-bit DT phi without its 4 lowest bits
  int dtPhi = wordPrecision ? ((digi.phi() >> 4) * 16) : digi.phi();
  stub.phiHw = angleConverter->getProcessorPhi(
      OMTFinputMaker::getProcessorPhiZero(config, iProcessor), procTyp, digi.scNum(), dtPhi);

  stub.etaHw = etaHw;

  if (iLayer == 0)
    stub.r = 431;  //round(431.175);  //MB1
  else if (iLayer == 2) {
    stub.r = 512;  //round(512.475);  //MB2
  } else if (iLayer == 4) {
    stub.r = 620;  //round(619.675);
    //MB3, it is different than in the phase-1, as in the phase-2 it is a middle of the DT chamber, not muon station
  }

  //in phase2, the phiB is 13 bits, and range is [-2, 2 rad] so 4 rad, 2^13 units/(4 rad) =  1^11/rad.
  //need to convert them to 512units==1rad (to use OLD PATTERNS...)
  stub.phiBHw = digi.phiBend() * config->dtPhiBUnitsRad() / 2048;
  //the cut if (stub.qualityHw >= config.getMinDtPhiBQuality()) is done in the ProcessorBase<GoldenPatternType>::restrictInput
  //as is is done like that in the firmware

  // need to shift 20-BX to roll-back the shift introduced by the DT TPs
  stub.bx = digi.bxNum() - config->dtBxShift();
  //stub.timing = digi.getTiming(); //TODO what about sub-bx timing, is is available?

  stub.logicLayer = iLayer;
  stub.detId = detid;

  OmtfName board(iProcessor, config);
  LogTrace("l1tOmtfEventPrint") << board.name() << " L1Phase2MuDTPhDigi: detid " << detid << " digi "
                                << " whNum " << digi.whNum() << " scNum " << digi.scNum() << " stNum " << digi.stNum()
                                << " slNum " << digi.slNum() << " quality " << digi.quality() << " rpcFlag "
                                << digi.rpcFlag() << " phi " << digi.phi() << " phiBend " << digi.phiBend() << " bx "
                                << digi.bxNum() << std::endl;
  LogTrace("l1tOmtfEventPrint") << board.name() << " stub: detid " << detid << " phi " << stub.phiHw << " eta "
                                << stub.etaHw << " phiB " << stub.phiBHw << " bx " << stub.bx << " quality "
                                << stub.qualityHw << " logicLayer " << stub.logicLayer << std::endl;
  OMTFinputMaker::addStub(config, muonStubsInLayers, iLayer, iInput, stub);
  if (outStub)
    *outStub = stub;
  return true;
}

void DtPhase2DigiToStubsConverterOmtf::addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers,
                                                     const L1Phase2MuDTThDigi& thetaDigi,
                                                     unsigned int iProcessor,
                                                     l1t::tftype procTyp) {
  //in the Phase1 omtf the theta stubs are merged with the phi in the addDTphiDigi
  //TODO implement if needed
}

bool DtPhase2DigiToStubsConverterOmtf::acceptDigi(const DTChamberId& dTChamberId,
                                                  unsigned int iProcessor,
                                                  l1t::tftype procType) {
  return OMTFinputMaker::acceptDtDigi(config, dTChamberId, iProcessor, procType);
}

InputMakerPhase2::InputMakerPhase2(const edm::ParameterSet& edmParameterSet,
                                   const MuStubsInputTokens& muStubsInputTokens,
                                   const MuStubsPhase2InputTokens& muStubsPhase2InputTokens,
                                   const OMTFConfiguration* config,
                                   std::unique_ptr<OmtfAngleConverter> angleConverter)
    : OMTFinputMaker(edmParameterSet, muStubsInputTokens, config, std::move(angleConverter)) {
  edm::LogImportant("OMTFReconstruction") << "constructing InputMakerPhase2" << std::endl;

  if (edmParameterSet.exists("usePhase2DTPrimitives") && edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
    if (edmParameterSet.getParameter<bool>("dropDTPrimitives") != true)
      throw cms::Exception(
          "L1TMuonOverlapPhase2 InputMakerPhase2::InputMakerPhase2 usePhase2DTPrimitives is true, but dropDTPrimitives "
          "is not true");
    //if the Phase2DTPrimitives are used, then the phase1 DT primitives should be dropped
    edm::LogImportant("OMTFReconstruction") << " using Phase2 DT trigger primitives" << std::endl;

    //dtPairPolicy: use the phi-theta pairs (DTTrigPhase2PairsProd) instead of separate phi and theta digis
    bool dtPairPolicy = edmParameterSet.exists("dtPairPolicy") && edmParameterSet.getParameter<bool>("dtPairPolicy");
    bool dtPairEtaWindowCheck =
        !edmParameterSet.exists("dtPairEtaWindowCheck") || edmParameterSet.getParameter<bool>("dtPairEtaWindowCheck");
    bool dtPairWordPrecision =
        edmParameterSet.exists("dtPairWordPrecision") && edmParameterSet.getParameter<bool>("dtPairWordPrecision");
    bool dtPairFixedEta = edmParameterSet.exists("dtPairFixedEta") && edmParameterSet.getParameter<bool>("dtPairFixedEta");
    if (dtPairPolicy)
      edm::LogImportant("OMTFReconstruction") << " dtPairPolicy: DT stubs and eta from phi-theta pairs" << std::endl;

    digiToStubsConverters.emplace_back(std::make_unique<DtPhase2DigiToStubsConverterOmtf>(
        config,
        dynamic_cast<OmtfPhase2AngleConverter*>(this->angleConverter.get()),
        muStubsPhase2InputTokens.inputTokenDtPh,
        muStubsPhase2InputTokens.inputTokenDtTh,
        muStubsPhase2InputTokens.inputTokenDtPairs,
        dtPairPolicy,
        dtPairEtaWindowCheck,
        dtPairFixedEta,
        dtPairWordPrecision));
  }
}
