/*
 * InputMakerPhase2.h
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#ifndef L1Trigger_L1TMuonOverlapPhase2_InputMakerPhase2_h
#define L1Trigger_L1TMuonOverlapPhase2_InputMakerPhase2_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTThContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTExtPhiThetaPairContainer.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfPhase2AngleConverter.h"

struct MuStubsPhase2InputTokens {
  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh;
  edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh;
  //only valid when dtPairPolicy is true
  edm::EDGetTokenT<L1Phase2MuDTExtPhiThetaPairContainer> inputTokenDtPairs;
};

class DtPhase2DigiToStubsConverter : public DigiToStubsConverterBase {
public:
  DtPhase2DigiToStubsConverter(const OMTFConfiguration* config,
                               edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh,
                               edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh,
                               edm::EDGetTokenT<L1Phase2MuDTExtPhiThetaPairContainer> inputTokenDtPairs = {},
                               bool dtPairPolicy = false)
      : config(config),
        dtPairPolicy(dtPairPolicy),
        inputTokenDtPh(inputTokenDtPh),
        inputTokenDtTh(inputTokenDtTh),
        inputTokenDtPairs(inputTokenDtPairs) {}

  ~DtPhase2DigiToStubsConverter() override {}

  void loadDigis(const edm::Event& event) override;

  void makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                 unsigned int iProcessor,
                 l1t::tftype procTyp,
                 int bxFrom,
                 int bxTo,
                 std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers,
                 XmlIOCache& xmlCache) override;

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                            const L1Phase2MuDTPhDigi& digi,
                            const L1Phase2MuDTThContainer* dtThDigis,
                            unsigned int iProcessor,
                            l1t::tftype procTyp) = 0;

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers,
                             const L1Phase2MuDTThDigi& thetaDigi,
                             unsigned int iProcessor,
                             l1t::tftype procTyp) = 0;

  //dtPairPolicy mode: one phi-theta pair per slot; the eta comes from the pair's own theta
  virtual void addDTpair(MuonStubPtrs2D& muonStubsInLayers,
                         const L1Phase2MuDTExtPhiThetaPair& pair,
                         unsigned int iProcessor,
                         l1t::tftype procTyp) {}

  virtual bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) {
    return true;
  }

protected:
  const OMTFConfiguration* config = nullptr;
  bool mergePhiAndTheta = true;
  bool dtPairPolicy = false;

  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh;
  edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh;
  edm::EDGetTokenT<L1Phase2MuDTExtPhiThetaPairContainer> inputTokenDtPairs;

  edm::Handle<L1Phase2MuDTPhContainer> dtPhDigis;
  edm::Handle<L1Phase2MuDTThContainer> dtThDigis;
  edm::Handle<L1Phase2MuDTExtPhiThetaPairContainer> dtPairs;

  int bunchCrossing = 0;
};

class DtPhase2DigiToStubsConverterOmtf : public DtPhase2DigiToStubsConverter {
public:
  DtPhase2DigiToStubsConverterOmtf(const OMTFConfiguration* config,
                                   const OmtfPhase2AngleConverter* angleConverter,
                                   edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh,
                                   edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh,
                                   edm::EDGetTokenT<L1Phase2MuDTExtPhiThetaPairContainer> inputTokenDtPairs = {},
                                   bool dtPairPolicy = false,
                                   bool dtPairEtaWindowCheck = true,
                                   bool dtPairFixedEta = false)
      : DtPhase2DigiToStubsConverter(config, inputTokenDtPh, inputTokenDtTh, inputTokenDtPairs, dtPairPolicy),
        angleConverter(angleConverter),
        dtPairEtaWindowCheck(dtPairEtaWindowCheck),
        dtPairFixedEta(dtPairFixedEta) {}

  ~DtPhase2DigiToStubsConverterOmtf() override = default;

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                    const L1Phase2MuDTPhDigi& digi,
                    const L1Phase2MuDTThContainer* dtThDigis,
                    unsigned int iProcessor,
                    l1t::tftype procTyp) override;

  void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers,
                     const L1Phase2MuDTThDigi& thetaDigi,
                     unsigned int iProcessor,
                     l1t::tftype procTyp) override;

  void addDTpair(MuonStubPtrs2D& muonStubsInLayers,
                 const L1Phase2MuDTExtPhiThetaPair& pair,
                 unsigned int iProcessor,
                 l1t::tftype procTyp) override;

  bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) override;

private:
  //common part of addDTphiDigi and addDTpair
  void addPhiStub(MuonStubPtrs2D& muonStubsInLayers,
                  const L1Phase2MuDTPhDigi& digi,
                  int etaHw,
                  unsigned int iProcessor,
                  l1t::tftype procTyp);

  const OmtfPhase2AngleConverter* angleConverter;
  bool dtPairEtaWindowCheck = true;
  bool dtPairFixedEta = false;
};

class InputMakerPhase2 : public OMTFinputMaker {
public:
  InputMakerPhase2(const edm::ParameterSet& edmParameterSet,
                   const MuStubsInputTokens& muStubsInputTokens,
                   const MuStubsPhase2InputTokens& muStubsPhase2InputTokens,
                   const OMTFConfiguration* config,
                   std::unique_ptr<OmtfAngleConverter> angleConverter);

  ~InputMakerPhase2() override = default;

  //the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers,
                            const L1Phase2MuDTPhDigi& digi,
                            const L1Phase2MuDTPhContainer* dtThDigis,
                            unsigned int iProcessor,
                            l1t::tftype procTyp) {}
};

#endif /* L1Trigger_L1TMuonOverlapPhase2_InputMakerPhase2_h */
