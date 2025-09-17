#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <iterator>
#include <utility>

// Helper function to map hwNumber to hwName based on the layer mapping
std::string getHwNameFromHwNumber(unsigned int hwNumber) {
  switch (hwNumber) {
    case 101: return "MB1";
    case 1101: return "MB1b";
    case 102: return "MB2";
    case 1102: return "MB2b";
    case 103: return "MB3";
    case 1103: return "MB3b";
    case 201: return "ME1/3";
    case 202: return "ME2/2";
    case 203: return "ME3/2";
    case 2011: return "ME1/2";
    case 301: return "RB1in";
    case 302: return "RB1out";
    case 303: return "RB2in";
    case 304: return "RB2out";
    case 305: return "RB3";
    case 311: return "RE1/3";
    case 312: return "RE2/3";
    case 313: return "RE3/3";
    default: return "Unknown";
  }
}

/////////////////////////////////////
void DtDigiToStubsConverter::loadDigis(const edm::Event& event) {
  event.getByToken(inputTokenDtPh, dtPhDigis);
  event.getByToken(inputTokenDtTh, dtThDigis);
}

void DtDigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                       unsigned int iProcessor,
                                       l1t::tftype procTyp,
                                       int bxFrom,
                                       int bxTo,
                                       std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  boost::property_tree::ptree procDataTree;
  
  for (const auto& digiIt : *dtPhDigis->getContainer()) {
    DTChamberId detid(digiIt.whNum(), digiIt.stNum(), digiIt.scNum() + 1);

    ///Check it the data fits into given processor input range
    if (!acceptDigi(detid, iProcessor, procTyp))
      continue;

    if (digiIt.bxNum() >= bxFrom && digiIt.bxNum() <= bxTo) {
      addDTphiDigi(muonStubsInLayers, digiIt, dtThDigis.product(), iProcessor, procTyp);
    }
  }

  if (!mergePhiAndTheta) {
    for (auto& thetaDigi : (*(dtThDigis->getContainer()))) {
      if (thetaDigi.bxNum() >= bxFrom && thetaDigi.bxNum() <= bxTo) {
        addDTetaStubs(muonStubsInLayers, thetaDigi, iProcessor, procTyp);
      }
    }
  }
  
  // Notify observers
  for (auto& observer : observers) {
    observer->addProcesorData("linkData", procDataTree);
  }
  //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" iProcessor "<<iProcessor<<std::endl;
}
///////////////////////////////////////
///////////////////////////////////////

void CscDigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                        unsigned int iProcessor,
                                        l1t::tftype procTyp,
                                        int bxFrom,
                                        int bxTo,
                                        std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  boost::property_tree::ptree procDataTree;

  auto chamber = cscDigis->begin();
  auto chend = cscDigis->end();
  for (; chamber != chend; ++chamber) {
    unsigned int rawid = (*chamber).first;
    ///Check it the data fits into given processor input range
    CSCDetId csc(rawid);
    if (!acceptDigi(csc, iProcessor, procTyp))
      continue;

    auto digi = (*chamber).second.first;
    auto dend = (*chamber).second.second;
    for (; digi != dend; ++digi) {
      ///Check if LCT trigger primitive has the right BX.
      int digiBx = digi->getBX() - config->cscLctCentralBx();

      if (digiBx >= bxFrom && digiBx <= bxTo) {
        addCSCstubs(muonStubsInLayers, rawid, *digi, iProcessor, procTyp);

        // === ADD XML DATA COLLECTION FOR CSV EXPORT ===
        auto& cscDigi = procDataTree.add_child("cscDigi", boost::property_tree::ptree());
        // Detector ID fields
        cscDigi.add("<xmlattr>.endcap", csc.endcap());
        cscDigi.add("<xmlattr>.station", csc.station());
        cscDigi.add("<xmlattr>.ring", csc.ring());
        cscDigi.add("<xmlattr>.chamber", csc.chamber());
        cscDigi.add("<xmlattr>.layer", csc.layer());
        
        // Core LCT data fields
        cscDigi.add("<xmlattr>.trknmb", digi->getTrknmb());
        cscDigi.add("<xmlattr>.valid", digi->isValid());
        cscDigi.add("<xmlattr>.quality", digi->getQuality());
        cscDigi.add("<xmlattr>.keywire", digi->getKeyWG());
        cscDigi.add("<xmlattr>.strip", digi->getStrip());
        cscDigi.add("<xmlattr>.pattern", digi->getPattern());
        cscDigi.add("<xmlattr>.bend", digi->getBend());
        cscDigi.add("<xmlattr>.bx", digi->getBX());
        cscDigi.add("<xmlattr>.mpclink", digi->getMPCLink());
        cscDigi.add("<xmlattr>.bx0", digi->getBX0());
        cscDigi.add("<xmlattr>.syncErr", digi->getSyncErr());
        cscDigi.add("<xmlattr>.cscID", digi->getCSCID());
        
        // Run-3 specific fields
        if (digi->isRun3()) {
          cscDigi.add("<xmlattr>.isRun3", true);
          cscDigi.add("<xmlattr>.quartStripBit", digi->getQuartStripBit());
          cscDigi.add("<xmlattr>.eighthStripBit", digi->getEighthStripBit());
          cscDigi.add("<xmlattr>.run3Pattern", digi->getRun3Pattern());
          cscDigi.add("<xmlattr>.slope", digi->getSlope());
          cscDigi.add("<xmlattr>.hmt", digi->getHMT());
        } else {
          cscDigi.add("<xmlattr>.isRun3", false);
        }
        
        // Additional computed fields
        cscDigi.add("<xmlattr>.fractionalStrip", digi->getFractionalStrip());
        cscDigi.add("<xmlattr>.fractionalSlope", digi->getFractionalSlope());
        cscDigi.add("<xmlattr>.clctPattern", digi->getCLCTPattern());
        cscDigi.add("<xmlattr>.stripType", digi->getStripType());
        cscDigi.add("<xmlattr>.bxData", digi->getBXData());
        cscDigi.add("<xmlattr>.type", digi->getType());
        
        // CSC angle conversion parameters
        CscConversionInfo convInfo = getCscConversionInfo(rawid, *digi, iProcessor, procTyp);
        cscDigi.add("<xmlattr>.offset", convInfo.offset);
        cscDigi.add("<xmlattr>.scale", convInfo.scale);
        cscDigi.add("<xmlattr>.order", convInfo.order);
        
        // === ADD GOLDEN STUB DATA ===
        // Check if a CSC stub was added to muonStubsInLayers
        const OMTFConfiguration* omtfConfig = dynamic_cast<const OMTFConfiguration*>(config);
        if (omtfConfig) {
          unsigned int hwNumber = omtfConfig->getLayerNumber(rawid);
          if (omtfConfig->getHwToLogicLayer().find(hwNumber) != omtfConfig->getHwToLogicLayer().end()) {
            unsigned int iLayer = omtfConfig->getHwToLogicLayer().at(hwNumber);
            if (iLayer < muonStubsInLayers.size()) {
              unsigned int iInput = OMTFinputMaker::getInputNumber(omtfConfig, rawid, iProcessor, procTyp);
              if (iInput < muonStubsInLayers[iLayer].size() && muonStubsInLayers[iLayer][iInput]) {
                auto& stub = muonStubsInLayers[iLayer][iInput];
                
                // Add golden stub data to XML
                auto& goldenStub = procDataTree.add_child("goldenStub", boost::property_tree::ptree());
                goldenStub.add("<xmlattr>.type", static_cast<int>(stub->type));
                goldenStub.add("<xmlattr>.logicLayer", stub->logicLayer);
                goldenStub.add("<xmlattr>.phiHw", stub->phiHw);
                goldenStub.add("<xmlattr>.etaHw", stub->etaHw);
                goldenStub.add("<xmlattr>.qualityHw", stub->qualityHw);
                goldenStub.add("<xmlattr>.phiBHw", stub->phiBHw);
                goldenStub.add("<xmlattr>.bx", stub->bx);
                goldenStub.add("<xmlattr>.timing", stub->timing);
                goldenStub.add("<xmlattr>.r", stub->r);
                goldenStub.add("<xmlattr>.detId", stub->detId);
                
                // Find hwNumber from logic layer and get hwName
                unsigned int stubHwNumber = omtfConfig->getLogicToHwLayer().at(stub->logicLayer);
                std::string hwName = getHwNameFromHwNumber(stubHwNumber);
                goldenStub.add("<xmlattr>.hwName", hwName);
              }
            }
          }
        }
      }
    }
  }

  // === NOTIFY OBSERVERS WITH XML DATA ===
  for (auto& observer : observers) {
    observer->addProcesorData("CSC", procDataTree);
  }
}

void RpcDigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                                        unsigned int iProcessor,
                                        l1t::tftype procTyp,
                                        int bxFrom,
                                        int bxTo,
                                        std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  //LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ <<" RPC HITS, processor : " << iProcessor<<" "<<std::endl;

  boost::property_tree::ptree procDataTree;

  const RPCDigiCollection& rpcDigiCollection = *rpcDigis;
  for (auto rollDigis : rpcDigiCollection) {
    RPCDetId roll = rollDigis.first;

    //debug
    //if(roll.region() != 0  &&  abs(roll.station()) >= 3 && roll.ring() == 1 )
    /*    {
      //iRPC
      for (auto pDigi=rollDigis.second.first; pDigi != rollDigis.second.second; pDigi++) {
        LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" roll "<<roll
            <<" strip "<<pDigi->strip()
            <<" hasX "<<pDigi->hasX()<<" coordinateX "<<pDigi->coordinateX()<<" hasY "<<pDigi->hasY()<<" coordinateY "<<pDigi->coordinateY()
            <<" bx "<<pDigi->bx()<<" time "<<pDigi->time()<<" irpc"<<std::endl;
      }
      //continue;
    }*/

    //LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ <<" roll "<<roll<<" "<<std::endl;

    if (!acceptDigi(roll, iProcessor, procTyp))
      continue;

    ///To find the clusters we have to copy the digis in chamber to sort them (not optimal).
    //  for (auto tdigi = rollDigis.second.first; tdigi != rollDigis.second.second; tdigi++) { std::cout << "RPC DIGIS: " << roll.rawId()<< " "<<roll<<" digi: " << tdigi->strip() <<" bx: " << tdigi->bx() << std::endl; }
    std::vector<RPCDigi> digisCopy;

    for (auto pDigi = rollDigis.second.first; pDigi != rollDigis.second.second; pDigi++) {
      if (pDigi->bx() >= bxFrom && pDigi->bx() <= bxTo) {
        digisCopy.push_back(*pDigi);

        // === ADD XML DATA COLLECTION FOR CSV EXPORT ===
        auto& rpcDigi = procDataTree.add_child("rpcDigi", boost::property_tree::ptree());
        // Detector ID fields
        rpcDigi.add("<xmlattr>.rpcID", roll.rawId());
        rpcDigi.add("<xmlattr>.region", roll.region());
        rpcDigi.add("<xmlattr>.ring", roll.ring());
        rpcDigi.add("<xmlattr>.station", roll.station());
        rpcDigi.add("<xmlattr>.sector", roll.sector());
        rpcDigi.add("<xmlattr>.layer", roll.layer());
        rpcDigi.add("<xmlattr>.subsector", roll.subsector());
        rpcDigi.add("<xmlattr>.roll", roll.roll());
        
        // Core RPC digi fields
        rpcDigi.add("<xmlattr>.strip", pDigi->strip());
        rpcDigi.add("<xmlattr>.bx", pDigi->bx());
        
        // Timing and spatial information
        rpcDigi.add("<xmlattr>.time", pDigi->time());
        rpcDigi.add("<xmlattr>.coordinateX", pDigi->coordinateX());
        rpcDigi.add("<xmlattr>.coordinateY", pDigi->coordinateY());
        rpcDigi.add("<xmlattr>.deltaTime", pDigi->deltaTime());
        rpcDigi.add("<xmlattr>.deltaX", pDigi->deltaX());
        rpcDigi.add("<xmlattr>.deltaY", pDigi->deltaY());
        
        // Flags for available information
        rpcDigi.add("<xmlattr>.hasTime", pDigi->hasTime());
        rpcDigi.add("<xmlattr>.hasX", pDigi->hasX());
        rpcDigi.add("<xmlattr>.hasY", pDigi->hasY());
        rpcDigi.add("<xmlattr>.isPseudoDigi", pDigi->isPseudoDigi());
      }
    }

    std::vector<RpcCluster> clusters = rpcClusterization->getClusters(roll, digisCopy);

    for (auto& cluster : clusters) {
      addRPCstub(muonStubsInLayers, roll, cluster, iProcessor, procTyp);
      
      // === ADD GOLDEN STUB DATA FOR RPC ===
      // Check if an RPC stub was added to muonStubsInLayers
      const OMTFConfiguration* omtfConfig = dynamic_cast<const OMTFConfiguration*>(config);
      if (omtfConfig) {
        unsigned int hwNumber = omtfConfig->getLayerNumber(roll.rawId());
        if (omtfConfig->getHwToLogicLayer().find(hwNumber) != omtfConfig->getHwToLogicLayer().end()) {
          unsigned int iLayer = omtfConfig->getHwToLogicLayer().at(hwNumber);
          if (iLayer < muonStubsInLayers.size()) {
            unsigned int iInput = OMTFinputMaker::getInputNumber(omtfConfig, roll.rawId(), iProcessor, procTyp);
            if (iInput < muonStubsInLayers[iLayer].size() && muonStubsInLayers[iLayer][iInput]) {
              auto& stub = muonStubsInLayers[iLayer][iInput];
              
              // Add golden stub data to XML
              auto& goldenStub = procDataTree.add_child("goldenStub", boost::property_tree::ptree());
              goldenStub.add("<xmlattr>.type", static_cast<int>(stub->type));
              goldenStub.add("<xmlattr>.logicLayer", stub->logicLayer);
              goldenStub.add("<xmlattr>.phiHw", stub->phiHw);
              goldenStub.add("<xmlattr>.etaHw", stub->etaHw);
              goldenStub.add("<xmlattr>.qualityHw", stub->qualityHw);
              goldenStub.add("<xmlattr>.phiBHw", stub->phiBHw);
              goldenStub.add("<xmlattr>.bx", stub->bx);
              goldenStub.add("<xmlattr>.timing", stub->timing);
              goldenStub.add("<xmlattr>.r", stub->r);
              goldenStub.add("<xmlattr>.detId", stub->detId);
              
              // Find hwNumber from logic layer and get hwName
              unsigned int stubHwNumber = omtfConfig->getLogicToHwLayer().at(stub->logicLayer);
              std::string hwName = getHwNameFromHwNumber(stubHwNumber);
              goldenStub.add("<xmlattr>.hwName", hwName);
            }
          }
        }
      }
    }
  }

  //removing the RPC stubs that were mark as dropped in the RpcDigiToStubsConverterOmtf::addRPCstub
  //10 is the first RPC layer
  for (unsigned int iLayer = 10; iLayer < muonStubsInLayers.size(); iLayer++) {
    for (unsigned int iInput = 0; iInput < muonStubsInLayers[iLayer].size(); iInput++) {
      if (muonStubsInLayers[iLayer][iInput] && muonStubsInLayers[iLayer][iInput]->type == MuonStub::RPC_DROPPED) {
        LogTrace("l1tOmtfEventPrint") << "RpcDigiToStubsConverter::makeStubs "
                                      << " iProcessor " << iProcessor << " procTyp " << procTyp
                                      << " dropping a stub iLayer " << iLayer << " iInput "
                                      << *(muonStubsInLayers[iLayer][iInput]) << std::endl;
        muonStubsInLayers[iLayer][iInput].reset();
      }
    }
  }

  // === NOTIFY OBSERVERS WITH XML DATA ===
  for (auto& observer : observers) {
    observer->addProcesorData("RPC", procDataTree);
  }
}

///////////////////////////////////////
///////////////////////////////////////
MuonStubMakerBase::MuonStubMakerBase(const ProcConfigurationBase* procConf) : config(procConf), rpcClusterization() {}

///////////////////////////////////////
///////////////////////////////////////
void MuonStubMakerBase::initialize(const edm::ParameterSet& edmCfg,
                                   const edm::EventSetup& es,
                                   const MuonGeometryTokens& muonGeometryTokens) {
  rpcClusterization.configure(
      config->getRpcMaxClusterSize(), config->getRpcMaxClusterCnt(), config->getRpcDropAllClustersIfMoreThanMax());
}
///////////////////////////////////////
///////////////////////////////////////
MuonStubMakerBase::~MuonStubMakerBase() {}
///////////////////////////////////////
///////////////////////////////////////

void MuonStubMakerBase::loadAndFilterDigis(const edm::Event& event) {
  for (auto& digiToStubsConverter : digiToStubsConverters)
    digiToStubsConverter->loadDigis(event);
}

void MuonStubMakerBase::buildInputForProcessor(MuonStubPtrs2D& muonStubsInLayers,
                                               unsigned int iProcessor,
                                               l1t::tftype procTyp,
                                               int bxFrom,
                                               int bxTo,
                                               std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  //LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ << " iProcessor " << iProcessor << " preocType "
  //                              << procTyp << std::endl;

  for (auto& digiToStubsConverter : digiToStubsConverters)
    digiToStubsConverter->makeStubs(muonStubsInLayers, iProcessor, procTyp, bxFrom, bxTo, observers);
}
