/*
 * MuonStubInputMaker.h
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef L1TkMuonBayes_MuonStubInputMaker_H_
#define L1TkMuonBayes_MuonStubInputMaker_H_

#include "L1Trigger/L1TkMuonBayes/interface/TkMuBayesProcConfig.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/AngleConverterBase.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"

namespace edm {
  class EventSetup;
}

class DtDigiToStubsConverterTkMu : public DtDigiToStubsConverter {
public:
  DtDigiToStubsConverterTkMu(const TkMuBayesProcConfig* config,
                             const AngleConverterBase* angleConverter,
                             edm::EDGetTokenT<L1MuDTChambPhContainer> inputTokenDtPh,
                             edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh)
      : DtDigiToStubsConverter(inputTokenDtPh, inputTokenDtTh), config(config), angleConverter(angleConverter){};

  ~DtDigiToStubsConverterTkMu() override{};

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                    const L1MuDTChambPhDigi& digi,
                    const L1MuDTChambThContainer* dtThDigis,
                    unsigned int iProcessor,
                    l1t::tftype procTyp) override;

  void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers,
                     const L1MuDTChambThDigi& thetaDigi,
                     unsigned int iProcessor,
                     l1t::tftype procTyp) override;

  static uint32_t getLayerNumber(const TkMuBayesProcConfig* config, const DTChamberId& detid, bool eta = false);
private:
  const TkMuBayesProcConfig* config = nullptr;
  const AngleConverterBase* angleConverter = nullptr;
};

class DtPhase2DigiToStubsConverterTkMu : public DtPhase2DigiToStubsConverter {
public:
  DtPhase2DigiToStubsConverterTkMu(const TkMuBayesProcConfig* config,
                                   const AngleConverterBase* angleConverter,
                                   edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh,
                                   edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh)
      : DtPhase2DigiToStubsConverter(inputTokenDtPh, inputTokenDtTh), config(config), angleConverter(angleConverter){};

  ~DtPhase2DigiToStubsConverterTkMu() override{};

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                    const L1Phase2MuDTPhDigi& digi,
                    const L1MuDTChambThContainer* dtThDigis,
                    unsigned int iProcessor,
                    l1t::tftype procTyp) override;

  //TODO no phase2 DTChambThDigi yet, so we are using the phase1
  void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers,
                     const L1MuDTChambThDigi& thetaDigi,
                     unsigned int iProcessor,
                     l1t::tftype procTyp) override;

private:
  const TkMuBayesProcConfig* config = nullptr;
  const AngleConverterBase* angleConverter;
};

class CscDigiToStubsConverterTkMu : public CscDigiToStubsConverter {
public:
  CscDigiToStubsConverterTkMu(const TkMuBayesProcConfig* config,
                              const AngleConverterBase* angleConverter,
                              edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCsc)
      : CscDigiToStubsConverter(config, inputTokenCsc), config(config), angleConverter(angleConverter){};

  ~CscDigiToStubsConverterTkMu() override{};

  //can add both phi and eta stubs
  void addCSCstubs(MuonStubPtrs2D& muonStubsInLayers,
                   unsigned int rawid,
                   const CSCCorrelatedLCTDigi& digi,
                   unsigned int iProcessor,
                   l1t::tftype procTyp) override;

  virtual uint32_t getLayerNumber(const CSCDetId& detid, bool eta = false) const;

private:
  const TkMuBayesProcConfig* config = nullptr;
  const AngleConverterBase* angleConverter;
};

class RpcDigiToStubsConverterTkMu : public RpcDigiToStubsConverter {
public:
  RpcDigiToStubsConverterTkMu(const TkMuBayesProcConfig* config,
                              const AngleConverterBase* angleConverter,
                              const RpcClusterization* rpcClusterization,
                              edm::EDGetTokenT<RPCDigiCollection> inputTokenRpc)
      : RpcDigiToStubsConverter(config, inputTokenRpc, rpcClusterization),
        config(config),
        angleConverter(angleConverter){};

  ~RpcDigiToStubsConverterTkMu() override{};

  void addRPCstub(MuonStubPtrs2D& muonStubsInLayers,
                  const RPCDetId& roll,
                  const RpcCluster& cluster,
                  unsigned int iProcessor,
                  l1t::tftype procTyp) override;

  virtual uint32_t getLayerNumber(const RPCDetId& detid) const;

private:
  const TkMuBayesProcConfig* config = nullptr;
  const AngleConverterBase* angleConverter;
};

class MuonStubInputMaker : public MuonStubMakerBase {
public:
  MuonStubInputMaker(const edm::ParameterSet& edmParameterSet,
                         MuStubsInputTokens& muStubsInputTokens,
                         edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2,
                         const TkMuBayesProcConfig* config,
                         AngleConverterBase* angleConv);

  ~MuonStubInputMaker() override;

  void initialize(const edm::ParameterSet& edmCfg, const edm::EventSetup& es) override;

  static void addStub(const TkMuBayesProcConfig* config,
                      MuonStubPtrs2D& muonStubsInLayers,
                      unsigned int iLayer,
                      MuonStub& stub);

private:
  const TkMuBayesProcConfig* config;

  std::unique_ptr<AngleConverterBase> angleConverter;
};

#endif /* L1TkMuonBayes_MuonStubInputMaker_H_ */
