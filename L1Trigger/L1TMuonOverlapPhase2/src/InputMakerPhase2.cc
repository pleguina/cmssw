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
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"

#include <iostream>
#include <algorithm>

/////////////////////////////////////
void DtPhase2DigiToStubsConverter::loadDigis(const edm::Event& event) {
  event.getByToken(inputTokenDtPh, dtPhDigis);
  event.getByToken(inputTokenDtTh, dtThDigis);
}

void DtPhase2DigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                             unsigned int iProcessor,
                                             l1t::tftype procTyp,
                                             int bxFrom,
                                             int bxTo,
                                             std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  boost::property_tree::ptree dtDigisTree;
  boost::property_tree::ptree dtStubsTree;
  
  // Maps to track DT digi ordering from the same chamber
  std::map<DTChamberId, int> chamberPhiDigiOrder;
  std::map<DTChamberId, int> chamberThDigiOrder;

  for (const auto& digiIt : *dtPhDigis->getContainer()) {
    DTChamberId detid(digiIt.whNum(), digiIt.stNum(), digiIt.scNum() + 1);

    ///Check it the data fits into given processor input range
    if (!acceptDigi(detid, iProcessor, procTyp))
      continue;

    // HACK for Phase-2  (DT TPs are centered in bX=20)
    if (digiIt.bxNum() - 20 >= bxFrom && digiIt.bxNum() - 20 <= bxTo) {
      addDTphiDigi(muonStubsInLayers, digiIt, dtThDigis.product(), iProcessor, procTyp);

      // Track ordering for multiple digis from the same chamber
      int chamberOrder = chamberPhiDigiOrder[detid]++;

      // Get hardware name for this DT chamber
      std::string hwName = getHwNameForDtChamber(detid);

      auto& dtP2PhiDigi = dtDigisTree.add_child("dtP2PhiDigi", boost::property_tree::ptree());
      dtP2PhiDigi.add("<xmlattr>.dtID", detid.rawId());
      dtP2PhiDigi.add("<xmlattr>.whNum", digiIt.whNum());
      dtP2PhiDigi.add("<xmlattr>.scNum", digiIt.scNum());
      dtP2PhiDigi.add("<xmlattr>.stNum", digiIt.stNum());
      dtP2PhiDigi.add("<xmlattr>.slNum", digiIt.slNum());
      dtP2PhiDigi.add("<xmlattr>.quality", digiIt.quality());
      dtP2PhiDigi.add("<xmlattr>.rpcFlag", digiIt.rpcFlag());
      dtP2PhiDigi.add("<xmlattr>.phi", digiIt.phi());
      dtP2PhiDigi.add("<xmlattr>.phiBend", digiIt.phiBend());
      dtP2PhiDigi.add("<xmlattr>.bx", digiIt.bxNum() - 20);  // ADD BX for CSV export
      dtP2PhiDigi.add("<xmlattr>.chamberOrder", chamberOrder);
      dtP2PhiDigi.add("<xmlattr>.hwName", hwName);  // ADD hwName attribute
    }
  }

  for (auto& thetaDigi : (*(dtThDigis->getContainer()))) {
    if (thetaDigi.bxNum() - 20 >= bxFrom && thetaDigi.bxNum() - 20 <= bxTo) {
      if (!mergePhiAndTheta) {
        addDTetaStubs(muonStubsInLayers, thetaDigi, iProcessor, procTyp);
      }

      // Track ordering for multiple digis from the same chamber
      DTChamberId thetaDetid(thetaDigi.whNum(), thetaDigi.stNum(), thetaDigi.scNum() + 1);
      int chamberOrder = chamberThDigiOrder[thetaDetid]++;

      // Get hardware name for this DT chamber
      std::string hwName = getHwNameForDtChamber(thetaDetid);

      auto& dtP2ThDigi = dtDigisTree.add_child("dtP2ThDigi", boost::property_tree::ptree());
      dtP2ThDigi.add("<xmlattr>.dtID", DTChamberId(thetaDigi.whNum(), thetaDigi.stNum(), thetaDigi.scNum() + 1).rawId());
      dtP2ThDigi.add("<xmlattr>.whNum", thetaDigi.whNum());
      dtP2ThDigi.add("<xmlattr>.scNum", thetaDigi.scNum());
      dtP2ThDigi.add("<xmlattr>.stNum", thetaDigi.stNum());
      dtP2ThDigi.add("<xmlattr>.quality", thetaDigi.quality());
      dtP2ThDigi.add("<xmlattr>.rpcFlag", thetaDigi.rpcFlag());
      dtP2ThDigi.add("<xmlattr>.k", thetaDigi.k());
      dtP2ThDigi.add("<xmlattr>.z", thetaDigi.z());
      dtP2ThDigi.add("<xmlattr>.bx", thetaDigi.bxNum() - 20);  // ADD BX for CSV export
      dtP2ThDigi.add("<xmlattr>.chamberOrder", chamberOrder);
      dtP2ThDigi.add("<xmlattr>.hwName", hwName);  // ADD hwName attribute
    }
  }

  // Capture DT golden stubs specifically (filter by stub type)
  // Track order per chamber (detId)
  std::map<uint32_t, int> chamberStubOrder;
  
  for (unsigned int iLayer = 0; iLayer < muonStubsInLayers.size(); iLayer++) {
    for (unsigned int iInput = 0; iInput < muonStubsInLayers[iLayer].size(); iInput++) {
      if (muonStubsInLayers[iLayer][iInput]) {
        const auto& stub = muonStubsInLayers[iLayer][iInput];
        // Only add DT stubs (DT types: DT_PHI, DT_THETA, DT_PHI_ETA, DT_HIT)
        if (stub->type == MuonStub::DT_PHI || stub->type == MuonStub::DT_THETA || 
            stub->type == MuonStub::DT_PHI_ETA || stub->type == MuonStub::DT_HIT) {
          
          // Track stub order per chamber
          uint32_t detId = stub->detId;
          int stubOrder = chamberStubOrder[detId]++;
          
          auto& dtStub = dtStubsTree.add_child("DTstub", boost::property_tree::ptree());
          dtStub.add("<xmlattr>.type", static_cast<int>(stub->type));
          dtStub.add("<xmlattr>.logicLayer", stub->logicLayer);
          dtStub.add("<xmlattr>.inputNumber", iInput);
          dtStub.add("<xmlattr>.phiHw", stub->phiHw);
          dtStub.add("<xmlattr>.etaHw", stub->etaHw);
          dtStub.add("<xmlattr>.qualityHw", stub->qualityHw);
          dtStub.add("<xmlattr>.phiBHw", stub->phiBHw);
          dtStub.add("<xmlattr>.bx", stub->bx);
          dtStub.add("<xmlattr>.timing", stub->timing);
          dtStub.add("<xmlattr>.r", stub->r);
          dtStub.add("<xmlattr>.detId", stub->detId);
          dtStub.add("<xmlattr>.order", stubOrder);  // ADD stub order
          
          // Add DT chamber identification fields
          DTChamberId dtId(stub->detId);
          dtStub.add("<xmlattr>.wheel", dtId.wheel());
          dtStub.add("<xmlattr>.station", dtId.station());
          dtStub.add("<xmlattr>.sector", dtId.sector());
          
          // Note: sector_wrapped is not added here because it requires config access
          // It will be handled in derived class implementations that have config access
          
          // Add hwName mapping - use virtual method to get the name
          std::string hwName = getHwNameForStub(stub->logicLayer);
          dtStub.add("<xmlattr>.hwName", hwName);
        }
      }
    }
  }

  for (auto& obs : observers) {
    obs->addProcesorData("DTdigis", dtDigisTree);
    obs->addProcesorData("DTstubs", dtStubsTree);
  }
}

//dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
void DtPhase2DigiToStubsConverterOmtf::addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                                                    const L1Phase2MuDTPhDigi& digi,
                                                    const L1Phase2MuDTThContainer* dtThDigis,
                                                    unsigned int iProcessor,
                                                    l1t::tftype procTyp) {
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

  if (stub.qualityHw < config.getMinDtPhiQuality())
    return;

  unsigned int hwNumber = config.getLayerNumber(detid.rawId());
  if (config.getHwToLogicLayer().find(hwNumber) == config.getHwToLogicLayer().end())
    return;

  auto iter = config.getHwToLogicLayer().find(hwNumber);
  unsigned int iLayer = iter->second;
  unsigned int iInput = OMTFinputMaker::getInputNumber(&config, detid.rawId(), iProcessor, procTyp);
  //MuonStub& stub = muonStubsInLayers[iLayer][iInput];

  stub.type = MuonStub::DT_PHI_ETA;

  stub.phiHw = angleConverter.getProcessorPhi(
      OMTFinputMaker::getProcessorPhiZero(&config, iProcessor), procTyp, digi.scNum(), digi.phi());

  stub.etaHw = angleConverter.getGlobalEta(detid, dtThDigis, digi.bxNum() - 20);

  if (iLayer == 0)
    stub.r = 431.175;  //MB1
  else if (iLayer == 2) {
    stub.r = 512.475;  //MB2
  } else if (iLayer == 4) {
    stub.r = 620;  //round(619.675);
    //MB3, it is different than in the phase-1, as in the phase-2 it is a middle of the DT chamber, not muon station
  }

  //in phase2, the phiB is 13 bits, and range is [-2, 2 rad] so 4 rad, 2^13 units/(4 rad) =  1^11/rad.
  //need to convert them to 512units==1rad (to use OLD PATTERNS...)
  stub.phiBHw = digi.phiBend() * config.dtPhiBUnitsRad() / 2048;
  //the cut if (stub.qualityHw >= config.getMinDtPhiBQuality()) is done in the ProcessorBase<GoldenPatternType>::restrictInput
  //as is is done like that in the firmware

  // need to shift 20-BX to roll-back the shift introduced by the DT TPs
  stub.bx = digi.bxNum() - 20;
  //stub.timing = digi.getTiming(); //TODO what about sub-bx timing, is is available?

  stub.logicLayer = iLayer;
  stub.detId = detid;

  OmtfName board(iProcessor, &config);
  LogTrace("l1tOmtfEventPrint") << board.name() << " L1Phase2MuDTPhDigi: detid " << detid << " digi "
                                << " whNum " << digi.whNum() << " scNum " << digi.scNum() << " stNum " << digi.stNum()
                                << " slNum " << digi.slNum() << " quality " << digi.quality() << " rpcFlag "
                                << digi.rpcFlag() << " phi " << digi.phi() << " phiBend " << digi.phiBend()
                                << std::endl;
  LogTrace("l1tOmtfEventPrint") << board.name() << " stub: detid " << detid << " phi " << stub.phiHw << " eta "
                                << stub.etaHw << " phiB " << stub.phiBHw << " bx " << stub.bx << " quality "
                                << stub.qualityHw << " logicLayer " << stub.logicLayer << std::endl;
  
  OMTFinputMaker::addStub(&config, muonStubsInLayers, iLayer, iInput, stub);
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
  return OMTFinputMaker::acceptDtDigi(&config, dTChamberId, iProcessor, procType);
}

std::string DtPhase2DigiToStubsConverterOmtf::getHwNameForStub(unsigned int logicLayer) {
  // Use the same mapping logic as Phase 1
  unsigned int stubHwNumber = config.getLogicToHwLayer().at(logicLayer);
  return getHwNameFromHwNumber(stubHwNumber);
}

InputMakerPhase2::InputMakerPhase2(const edm::ParameterSet& edmParameterSet,
                                   MuStubsInputTokens& muStubsInputTokens,
                                   MuStubsPhase2InputTokens& muStubsPhase2InputTokens,
                                   const OMTFConfiguration* config,
                                   std::unique_ptr<OmtfAngleConverter> angleConverter)
    : OMTFinputMaker(edmParameterSet, muStubsInputTokens, config, std::move(angleConverter)) {
  edm::LogImportant("OMTFReconstruction") << "constructing InputMakerPhase2" << std::endl;

  // Phase2 Fix: Check dropRPCPrimitives parameter
  if (edmParameterSet.getParameter<bool>("dropRPCPrimitives")) {
    // Remove RPC converter if dropRPCPrimitives is true
    digiToStubsConverters.erase(
        std::remove_if(digiToStubsConverters.begin(), digiToStubsConverters.end(),
                       [](const std::unique_ptr<DigiToStubsConverterBase>& converter) {
                         return dynamic_cast<RpcDigiToStubsConverterOmtf*>(converter.get()) != nullptr;
                       }),
        digiToStubsConverters.end());
    edm::LogImportant("OMTFReconstruction") << " dropping RPC primitives in Phase2" << std::endl;
  }

  // Phase2 Fix: Check dropCSCPrimitives parameter
  if (edmParameterSet.getParameter<bool>("dropCSCPrimitives")) {
    // Remove CSC converter if dropCSCPrimitives is true
    digiToStubsConverters.erase(
        std::remove_if(digiToStubsConverters.begin(), digiToStubsConverters.end(),
                       [](const std::unique_ptr<DigiToStubsConverterBase>& converter) {
                         return dynamic_cast<CscDigiToStubsConverterOmtf*>(converter.get()) != nullptr;
                       }),
        digiToStubsConverters.end());
    edm::LogImportant("OMTFReconstruction") << " dropping CSC primitives in Phase2" << std::endl;
  }

  if (edmParameterSet.exists("usePhase2DTPrimitives") && edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
    if (edmParameterSet.getParameter<bool>("dropDTPrimitives") != true)
      throw cms::Exception(
          "L1TMuonOverlapPhase2 InputMakerPhase2::InputMakerPhase2 usePhase2DTPrimitives is true, but dropDTPrimitives "
          "is not true");
    //if the Phase2DTPrimitives are used, then the phase1 DT primitives should be dropped
    edm::LogImportant("OMTFReconstruction") << " using Phase2 DT trigger primitives" << std::endl;

    digiToStubsConverters.emplace_back(std::make_unique<DtPhase2DigiToStubsConverterOmtf>(
        config,
        dynamic_cast<OmtfPhase2AngleConverter*>(this->angleConverter.get()),
        muStubsPhase2InputTokens.inputTokenDtPh,
        muStubsPhase2InputTokens.inputTokenDtTh));
  }
}

std::string DtPhase2DigiToStubsConverterOmtf::getHwNameForDtChamber(const DTChamberId& detid) {
  unsigned int hwNumber = config.getLayerNumber(detid.rawId());
  return getHwNameFromHwNumber(hwNumber);
}

// Override makeStubs to add reference stub functionality
void DtPhase2DigiToStubsConverterOmtf::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                                 unsigned int iProcessor,
                                                 l1t::tftype procTyp,
                                                 int bxFrom,
                                                 int bxTo,
                                                 std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  // Call the base class implementation first
  DtPhase2DigiToStubsConverter::makeStubs(muonStubsInLayers, iProcessor, procTyp, bxFrom, bxTo, observers);
  
  // Now handle reference stubs specifically for OMTF - use global tree
  
  // Capture DT reference stubs specifically (filter by stub type)
  for (unsigned int iLayer = 0; iLayer < muonStubsInLayers.size(); iLayer++) {
    for (unsigned int iInput = 0; iInput < muonStubsInLayers[iLayer].size(); iInput++) {
      if (muonStubsInLayers[iLayer][iInput]) {
        const auto& stub = muonStubsInLayers[iLayer][iInput];
        // Only add DT stubs (DT types: DT_PHI, DT_THETA, DT_PHI_ETA, DT_HIT)
        if (stub->type == MuonStub::DT_PHI || stub->type == MuonStub::DT_THETA || 
            stub->type == MuonStub::DT_PHI_ETA || stub->type == MuonStub::DT_HIT) {
          
          // Check if this is a reference layer stub and add to reference collection
          const auto& refToLogicNumbers = config.getRefToLogicNumber();
          for (unsigned int iRefLayer = 0; iRefLayer < refToLogicNumbers.size(); ++iRefLayer) {
            if (refToLogicNumbers[iRefLayer] == (int)stub->logicLayer) {
              // This is a reference layer stub - add to reference collection
              DTChamberId dtId(stub->detId);
              
              // Calculate sector_wrapped based on sector number
              // Configuration: Proc 0: sectors {1,2,3,4,5}, Proc 1: {5,6,7,8,9}, Proc 2: {9,10,11,12,1}
              // Sectors 1, 5, and 9 are shared (overlap) between adjacent processors
              // Hardware expects indices 0-4 for each processor
              // Formula: sector_wrapped = sector - barrelMin[iProcessor]
              int sector = dtId.sector();
              int sector_wrapped = 0;
              
              // Get barrel minimum for this processor
              int aMin = config.getBarrelMin()[iProcessor];
              
              // Handle wrap-around for last processor (sectors 1,2 belong to proc 2)
              int aSector = sector;
              if (iProcessor == (config.nProcessors() - 1) && aSector < 3) {
                aSector += 12;  // 12 sectors total in barrel
              }
              
              // Calculate sector_wrapped using the standard formula
              sector_wrapped = aSector - aMin;
              
              // Add hwName mapping
              std::string hwName = getHwNameForStub(stub->logicLayer);
              
              // Add region classification based on phi value
              unsigned int logicRegion = calculateLogicRegion(stub->phiHw, iRefLayer, iInput, &config);
              
              // Use static method to add reference stub
              MuonStubMakerBase::addGlobalReferenceStub("DT", iProcessor, iRefLayer, stub->logicLayer, stub->phiHw, stub->phiBHw, stub->etaHw, 
                                                       stub->qualityHw, stub->detId, hwName, dtId.wheel(), dtId.station(), 
                                                       -1, -1, -1,  // CSC fields: not applicable for DT
                                                       dtId.sector(), sector_wrapped,  // DT fields: sector=dtId.sector(), sector_wrapped=sector_wrapped
                                                       logicRegion);
              break; // Found the reference layer, no need to continue
            }
          }
        }
      }
    }
  }

  // Reference stubs will be output by MuonStubMakerBase::flushReferenceStubs()
}


