/*
 * OmtfProcessorLut2d.h
 *
 *  Created on: May 14, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF2_OMTFPROCESSORLUT2D_H_
#define INTERFACE_OMTF2_OMTFPROCESSORLUT2D_H_

#include "L1Trigger/L1TMuonBayes/interface/MuonStubsInput.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/AlgoMuon2.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/Lut2d.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/PtModuleLut2D.h"

#include "L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

class IStubSelector {
public:
  virtual ~IStubSelector() {}

  //makes from the muonStubs a vector of MuonStubsInput if, but each MuonStubsInput in the vector contain at most one stub per layer
  virtual std::vector<MuonStubsInput> selelctStubs(const MuonStubsInput& muonStubs) = 0;
};

class SimpleStubSelector: public IStubSelector {
public:
  SimpleStubSelector(const ProcConfigurationBase* config): IStubSelector(), config(config) {}

  virtual ~SimpleStubSelector() {}

  virtual std::vector<MuonStubsInput> selelctStubs(const MuonStubsInput& muonStubs);
private:
  const ProcConfigurationBase* config;
};

class OmtfProcessorLut2d {
public:
  OmtfProcessorLut2d(OMTFConfiguration* omtfConfig);

  OmtfProcessorLut2d(OMTFConfiguration* omtfConfig, const std::shared_ptr<PtModuleLut2D>& ptModule);

  virtual ~OmtfProcessorLut2d();

  //virtual MuonStubsInput selectStubs(const MuonStubsInput& muonStubs);

  virtual AlgoMuon2s processStubs(const MuonStubsInput& muonStubs);

  //convert algo muon to outgoing Candidates
  virtual std::vector<l1t::RegionalMuonCand> getFinalcandidates(
                 unsigned int iProcessor, l1t::tftype mtfType, const AlgoMuon2s& algoCands);

  std::shared_ptr<PtModuleLut2D>& getPtModule() {
    return ptModule;
  }
private:
  OMTFConfiguration* config;

  std::unique_ptr<IStubSelector> stubSelector;

  std::shared_ptr<PtModuleLut2D> ptModule;

};

#endif /* INTERFACE_OMTF2_OMTFPROCESSORLUT2D_H_ */
