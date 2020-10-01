/*
 * MuonCorrelatorInputMaker.h
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_
#define MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_


#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "L1Trigger/L1TkMuonBayes/interface/MuCorrelatorConfig.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/AngleConverterBase.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"

namespace edm {
  class EventSetup;
}

class DtDigiToStubsConverterTkMu: public DtDigiToStubsConverter {
public:
  DtDigiToStubsConverterTkMu(const MuCorrelatorConfig* config, const AngleConverterBase* angleConverter, edm::EDGetTokenT<L1MuDTChambPhContainer> inputTokenDtPh, edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh):
    DtDigiToStubsConverter(inputTokenDtPh, inputTokenDtTh), config(config), angleConverter(angleConverter) {};

  virtual ~DtDigiToStubsConverterTkMu() {};

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambPhDigi& digi, const L1MuDTChambThContainer *dtThDigis,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
      unsigned int iProcessor, l1t::tftype procTyp);

private:
  const MuCorrelatorConfig* config =  nullptr;
  const AngleConverterBase* angleConverter =  nullptr;
};

class DtPhase2DigiToStubsConverterTkMu: public DtPhase2DigiToStubsConverter {
public:
  DtPhase2DigiToStubsConverterTkMu(const MuCorrelatorConfig* config, const AngleConverterBase* angleConverter, edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh, edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh):
    DtPhase2DigiToStubsConverter(inputTokenDtPh, inputTokenDtTh), config(config), angleConverter(angleConverter) {};

  virtual ~DtPhase2DigiToStubsConverterTkMu() {};

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1Phase2MuDTPhDigi& digi, const L1MuDTChambThContainer *dtThDigis,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
      unsigned int iProcessor, l1t::tftype procTyp);

private:
  const MuCorrelatorConfig* config =  nullptr;
  const AngleConverterBase* angleConverter;
};

class CscDigiToStubsConverterTkMu: public CscDigiToStubsConverter {
public:
  CscDigiToStubsConverterTkMu(const MuCorrelatorConfig* config, const AngleConverterBase* angleConverter, edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCsc):
    CscDigiToStubsConverter(config, inputTokenCsc), config(config), angleConverter(angleConverter) {};

  virtual ~CscDigiToStubsConverterTkMu() {};

  //can add both phi and eta stubs
  virtual void addCSCstubs(MuonStubPtrs2D& muonStubsInLayers, unsigned int rawid, const CSCCorrelatedLCTDigi& digi,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual uint32_t getLayerNumber(const CSCDetId& detid, bool eta = false) const;
private:
  const MuCorrelatorConfig* config =  nullptr;
  const AngleConverterBase* angleConverter;
};


class RpcDigiToStubsConverterTkMu: public RpcDigiToStubsConverter {
public:
  RpcDigiToStubsConverterTkMu(const MuCorrelatorConfig* config, const AngleConverterBase* angleConverter, const RpcClusterization* rpcClusterization, edm::EDGetTokenT<RPCDigiCollection> inputTokenRpc):
    RpcDigiToStubsConverter(config, inputTokenRpc, rpcClusterization), config(config), angleConverter(angleConverter) {};

  virtual ~RpcDigiToStubsConverterTkMu() {};

  virtual void addRPCstub(MuonStubPtrs2D& muonStubsInLayers, const RPCDetId& roll, const RpcCluster& cluster,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual uint32_t getLayerNumber(const RPCDetId& detid) const;
private:
  const MuCorrelatorConfig* config =  nullptr;
  const AngleConverterBase* angleConverter;
};


class MuCorrelatorInputMaker : public MuonStubMakerBase {
public:
  MuCorrelatorInputMaker(const edm::ParameterSet& edmParameterSet, MuStubsInputTokens& muStubsInputTokens,
      edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2, const MuCorrelatorConfig* config);

  virtual ~MuCorrelatorInputMaker();

  static void addStub(const MuCorrelatorConfig* config, MuonStubPtrs2D& muonStubsInLayers, unsigned int iLayer, MuonStub& stub);

  //logic layer numbers
  static uint32_t getLayerNumber(const MuCorrelatorConfig* config, const DTChamberId& detid, bool eta = false);
  //uint32_t getLayerNumber(const CSCDetId& detid, bool eta = false) const;
  //uint32_t getLayerNumber(const RPCDetId& detid) const;uint32_t getLayerNumber(const RPCDetId& detid) const;

private:

/*  virtual bool acceptDigi(uint32_t rawId, unsigned int iProcessor, l1t::tftype procType) {
    return true;
  }

  //the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambPhDigi& digi,
     const L1MuDTChambThContainer *dtThDigis,
     unsigned int iProcessor,
     l1t::tftype procTyp);

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
       unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addCSCstubs(MuonStubPtrs2D& muonStubsInLayers,  unsigned int rawid, const CSCCorrelatedLCTDigi& digi,
     unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addRPCstub(MuonStubPtrs2D& muonStubsInLayers, const RPCDetId& roll, const RpcCluster& cluster,
     unsigned int iProcessor, l1t::tftype procTyp);*/

  const MuCorrelatorConfig* config;

  AngleConverterBase angleConverter;
};

#endif /* INTERFACE_MUCORRELATOR_MUONCORRELATORINPUTMAKER_H_ */
