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
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/XmlIOCache.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfPhase2AngleConverter.h"

struct MuStubsPhase2InputTokens {
  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh;
  edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh;
};

class DtPhase2DigiToStubsConverter : public DigiToStubsConverterBase {
public:
  DtPhase2DigiToStubsConverter(edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh,
                               edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh)
      : inputTokenDtPh(inputTokenDtPh), inputTokenDtTh(inputTokenDtTh) {}

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

  virtual bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) {
    return true;
  }

  // Virtual method to get hwName for a DT chamber - implemented in derived class
  virtual std::string getHwNameForDtChamber(const DTChamberId& detid) { return ""; }

  // Virtual method to get hardware name from logic layer - to be implemented by derived classes
  virtual std::string getHwNameForStub(unsigned int logicLayer) {
    return "Unknown";  // Default implementation
  }

  // Virtual method to calculate sector_wrapped for DT - to be implemented by derived classes with config access
  virtual int calculateDtSectorWrapped(int sector, unsigned int iProcessor) {
    return -1;  // Default: no calculation (derived class should override)
  }

  // Virtual method to check if a layer is a bending layer - to be implemented by derived classes
  virtual bool isBendingLayer(unsigned int iLayer) {
    return false;  // Default: no bending layers (derived class should override)
  }

  // Virtual method to get minimum DT phiB quality for bending layers
  virtual int getMinDtPhiBQuality() {
    return 0;  // Default: no quality cut (derived class should override)
  }

protected:
  bool mergePhiAndTheta = true;

  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh;
  edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh;

  edm::Handle<L1Phase2MuDTPhContainer> dtPhDigis;
  edm::Handle<L1Phase2MuDTThContainer> dtThDigis;
};

class DtPhase2DigiToStubsConverterOmtf : public DtPhase2DigiToStubsConverter {
public:
  DtPhase2DigiToStubsConverterOmtf(const OMTFConfiguration* config,
                                   const OmtfPhase2AngleConverter* angleConverter,
                                   edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh,
                                   edm::EDGetTokenT<L1Phase2MuDTThContainer> inputTokenDtTh)
      : DtPhase2DigiToStubsConverter(inputTokenDtPh, inputTokenDtTh),
        config(*config),
        angleConverter(*angleConverter) {}

  ~DtPhase2DigiToStubsConverterOmtf() override = default;

  // Override makeStubs to add reference stub functionality
  void makeStubs(MuonStubPtrs2D& muonStubsInLayers,
                 unsigned int iProcessor,
                 l1t::tftype procTyp,
                 int bxFrom,
                 int bxTo,
                 std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers,
                 XmlIOCache& xmlCache) override;

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

  bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) override;

  // Implementation of hwName mapping for OMTF
  std::string getHwNameForStub(unsigned int logicLayer) override;
  
  // Implementation of hwName for DT chamber
  std::string getHwNameForDtChamber(const DTChamberId& detid) override;
  
  // Override to calculate sector_wrapped for DT using OMTF config
  int calculateDtSectorWrapped(int sector, unsigned int iProcessor) override;

  // Override to check if layer is a bending layer using OMTF config
  bool isBendingLayer(unsigned int iLayer) override;

  // Override to get minimum DT phiB quality from OMTF config
  int getMinDtPhiBQuality() override;

private:
  const OMTFConfiguration& config;
  const OmtfPhase2AngleConverter& angleConverter;
};

class InputMakerPhase2 : public OMTFinputMaker {
public:
  InputMakerPhase2(const edm::ParameterSet& edmParameterSet,
                   MuStubsInputTokens& muStubsInputTokens,
                   MuStubsPhase2InputTokens& muStubsPhase2InputTokens,
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
