/*
 * InputMakerPhase2.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"
#include <iostream>


/////////////////////////////////////
void DtPhase2DigiToStubsConverter::loadDigis(const edm::Event& event) {
  event.getByToken(inputTokenDtPh, dtPhDigis);
  event.getByToken(inputTokenDtTh, dtThDigis);
}

void DtPhase2DigiToStubsConverter::makeStubs(MuonStubPtrs2D& muonStubsInLayers, unsigned int iProcessor, l1t::tftype procTyp, int bxFrom, int bxTo) {
  if(!dtPhDigis)
    return;

  for (const auto& digiIt: *dtPhDigis->getContainer()) {
    DTChamberId detid(digiIt.whNum(),digiIt.stNum(),digiIt.scNum()+1);

    ///Check it the data fits into given processor input range
    if(!acceptDigi(detid, iProcessor, procTyp))
      continue;

    if (digiIt.bxNum() >= bxFrom && digiIt.bxNum() <= bxTo )
      addDTphiDigi(muonStubsInLayers, digiIt, dtThDigis.product() , iProcessor, procTyp);
  }

  if(!mergePhiAndTheta) {
    for(auto& thetaDigi: (*(dtThDigis->getContainer()) ) ) {
      if(thetaDigi.bxNum() >= bxFrom && thetaDigi.bxNum() <= bxTo) {
        addDTetaStubs(muonStubsInLayers, thetaDigi, iProcessor, procTyp);
      }
    }
  }
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iProcessor "<<iProcessor<<std::endl;
  //angleConverter->AngleConverterBase::getGlobalEta(dtThDigis, 0, 0);std::endl;
  //angleConverter->AngleConverterBase::getGlobalEta(dtThDigis, 0, 0);
}

//dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
void DtPhase2DigiToStubsConverterOmtf::addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1Phase2MuDTPhDigi& digi, const L1MuDTChambThContainer* dtThDigis,
    unsigned int iProcessor, l1t::tftype procTyp)
{
  DTChamberId detid(digi.whNum(), digi.stNum(), digi.scNum()+1);

  //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" OMTFinputMaker "<<" detid "<<detid<<endl;

// TODO FIXME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!
  if(digi.quality() ==  7 || digi.quality() < config->getMinDtPhiQuality()) //7 is empty digi, TODO update if the definition of the quality is changed
    return;

  unsigned int hwNumber = config->getLayerNumber(detid.rawId());
  if(config->getHwToLogicLayer().find(hwNumber) == config->getHwToLogicLayer().end())
    return;

  auto iter = config->getHwToLogicLayer().find(hwNumber);
  unsigned int iLayer = iter->second;
  unsigned int iInput = OMTFinputMaker::getInputNumber(config, detid.rawId(), iProcessor, procTyp);
  //MuonStub& stub = muonStubsInLayers[iLayer][iInput];
  MuonStub stub;

  stub.type = MuonStub::DT_PHI_ETA;
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iProcessor "<<iProcessor<<std::endl;
  stub.phiHw  =  angleConverter->getProcessorPhi(OMTFinputMaker::getProcessorPhiZero(config, iProcessor), procTyp, digi.scNum(), digi.phi());
  //stub.etaHw  =  angleConverter->getGlobalEta(digi, dtThDigis);
  stub.etaHw  =  angleConverter->getGlobalEta(detid, dtThDigis, digi.bxNum() );
  stub.phiBHw = digi.phiBend();
  stub.qualityHw = digi.quality();

  stub.bx = digi.bxNum(); //TODO sholdn't  it be BxCnt()?
  //stub.timing = digi.getTiming(); //TODO what about sub-bx timing, is is available?

  //stub.etaType = ?? TODO
  stub.logicLayer = iLayer;
  stub.detId = detid;

  OMTFinputMaker::addStub(config, muonStubsInLayers, iLayer, iInput, stub);
}

void DtPhase2DigiToStubsConverterOmtf::addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
    unsigned int iProcessor, l1t::tftype procTyp) {
  //in the Phase1 omtf the theta stubs are merged with the phi in the addDTphiDigi
  //TODO implement if needed
}

bool DtPhase2DigiToStubsConverterOmtf::acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) {
  return OMTFinputMaker::acceptDtDigi(config, dTChamberId, iProcessor, procType);
}

InputMakerPhase2::InputMakerPhase2(const edm::ParameterSet& edmParameterSet, MuStubsInputTokens& muStubsInputTokens,
    edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2, const OMTFConfiguration* config)
    :OMTFinputMaker(edmParameterSet, muStubsInputTokens, config)
{
	edm::LogImportant("OMTFReconstruction") << "constructing InputMakerPhase2" << std::endl;

/*  if(!edmParameterSet.getParameter<bool>("dropDTPrimitives"))
    digiToStubsConverters.emplace_back(std::make_unique<DtDigiToStubsConverterOmtf>(config, &angleConverter, muStubsInputTokens.inputTokenDtPh, muStubsInputTokens.inputTokenDtTh));

  if(!edmParameterSet.getParameter<bool>("dropCSCPrimitives"))
    digiToStubsConverters.emplace_back(std::make_unique<CscDigiToStubsConverterOmtf>(config, &angleConverter, muStubsInputTokens.inputTokenCSC));

  if(!edmParameterSet.getParameter<bool>("dropRPCPrimitives"))
    digiToStubsConverters.emplace_back(std::make_unique<RpcDigiToStubsConverterOmtf>(config, &angleConverter, &rpcClusterization, muStubsInputTokens.inputTokenRPC));*/

	if(edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
	  if(edmParameterSet.getParameter<bool>("dropDTPrimitives") != true)
	    throw cms::Exception("L1TMuonOverlapPhase2 InputMakerPhase2::InputMakerPhase2 usePhase2DTPrimitives is true, but dropDTPrimitives is not true");
	    //if the Phase2DTPrimitives are used, then the phase1 DT primitives should be dropped
	  digiToStubsConverters.emplace_back(std::make_unique<DtPhase2DigiToStubsConverterOmtf>(config, &angleConverter, inputTokenDTPhPhase2, muStubsInputTokens.inputTokenDtTh));
	}
}

InputMakerPhase2::~InputMakerPhase2() {
	// TODO Auto-generated destructor stub
}

