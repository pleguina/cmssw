/*
 * InputMakerPhase2.h
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_INPUTMAKERPHASE2_H_
#define INTERFACE_INPUTMAKERPHASE2_H_

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"


class DtPhase2DigiToStubsConverter: public DigiToStubsConverterBase {
public:
  DtPhase2DigiToStubsConverter(edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh, edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh):
    inputTokenDtPh(inputTokenDtPh), inputTokenDtTh(inputTokenDtTh) {};

  virtual ~DtPhase2DigiToStubsConverter() {};

  //virtual void initialize(const edm::ParameterSet& edmCfg, const edm::EventSetup& es, const ProcConfigurationBase* procConf) {};

  virtual void loadDigis(const edm::Event& event);

  virtual void makeStubs(MuonStubPtrs2D& muonStubsInLayers, unsigned int iProcessor, l1t::tftype procTyp, int bxFrom, int bxTo);

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1Phase2MuDTPhDigi& digi, const L1MuDTChambThContainer *dtThDigis,
      unsigned int iProcessor, l1t::tftype procTyp) = 0;

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
      unsigned int iProcessor, l1t::tftype procTyp) = 0;

  virtual bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType) {return true;}
protected:
  bool mergePhiAndTheta = true;

  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh;

  edm::Handle<L1Phase2MuDTPhContainer> dtPhDigis;
  edm::Handle<L1MuDTChambThContainer> dtThDigis;
};


class DtPhase2DigiToStubsConverterOmtf: public DtPhase2DigiToStubsConverter {
public:
  DtPhase2DigiToStubsConverterOmtf(const OMTFConfiguration* config, const OmtfAngleConverter* angleConverter, edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDtPh, edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDtTh):
    DtPhase2DigiToStubsConverter(inputTokenDtPh, inputTokenDtTh), config(config), angleConverter(angleConverter) {};

  virtual ~DtPhase2DigiToStubsConverterOmtf() {};

  //dtThDigis is provided as argument, because in the OMTF implementation the phi and eta digis are merged (even thought it is artificial)
  virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1Phase2MuDTPhDigi& digi, const L1MuDTChambThContainer *dtThDigis,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual void addDTetaStubs(MuonStubPtrs2D& muonStubsInLayers, const L1MuDTChambThDigi& thetaDigi,
      unsigned int iProcessor, l1t::tftype procTyp);

  virtual bool acceptDigi(const DTChamberId& dTChamberId, unsigned int iProcessor, l1t::tftype procType);

private:
  const OMTFConfiguration* config =  nullptr;
  const OmtfAngleConverter* angleConverter;
};


class InputMakerPhase2: public OMTFinputMaker {
public:
	InputMakerPhase2(const edm::ParameterSet& edmParameterSet, MuStubsInputTokens& muStubsInputTokens,
	    edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2, const OMTFConfiguration* config );

	virtual ~InputMakerPhase2();

	 //the phi and eta digis are merged (even thought it is artificial)
	 virtual void addDTphiDigi(MuonStubPtrs2D& muonStubsInLayers, const L1Phase2MuDTPhDigi& digi,
	    const L1Phase2MuDTPhContainer *dtThDigis,
	    unsigned int iProcessor,
	    l1t::tftype procTyp)  {}

private:
};

#endif /* INTERFACE_INPUTMAKERPHASE2_H_ */
