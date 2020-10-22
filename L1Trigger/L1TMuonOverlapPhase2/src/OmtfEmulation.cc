/*
 * OmtfEmulation.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfEmulation.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"
//#include "L1Trigger/L1TMuonOverlapPhase2/interface/PtAssignmentNN.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

OmtfEmulation::OmtfEmulation(const edm::ParameterSet& edmParameterSet,
                             MuStubsInputTokens& muStubsInputTokens,
                             edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2)
    : OMTFReconstruction(edmParameterSet, muStubsInputTokens), inputTokenDTPhPhase2(inputTokenDTPhPhase2) {}

OmtfEmulation::~OmtfEmulation() {}

void OmtfEmulation::beginJob() {
  if (edmParameterSet.exists("usePhase2DTPrimitives") && edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
    inputMaker.reset(new InputMakerPhase2(
        edmParameterSet, muStubsInputTokens, inputTokenDTPhPhase2, omtfConfig.get(), new OmtfPhase2AngleConverter()));
  } else {
    inputMaker = std::make_unique<OMTFinputMaker>(
        edmParameterSet, muStubsInputTokens, omtfConfig.get(), new OmtfAngleConverter());
  }
}

void OmtfEmulation::addObservers() {
  OMTFReconstruction::addObservers();

  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>*>(omtfProc.get());
  if (omtfProcGoldenPat) {
    /*if (edmParameterSet.exists("neuralNetworkFile")) {
      edm::LogImportant("OMTFReconstruction") << "constructing PtAssignmentNN" << std::endl;
      std::string neuralNetworkFile = edmParameterSet.getParameter<edm::FileInPath>("neuralNetworkFile").fullPath();
      omtfProcGoldenPat->setPtAssignment(new PtAssignmentNN(
          edmParameterSet, omtfConfig.get(), neuralNetworkFile));  //TODO change to dynamic_cast and check the type
    }*/

    /*    if(edmParameterSet.exists("patternsPtAssignment") && edmParameterSet.getParameter<bool>("patternsPtAssignment")) {
      //std::string rootFileName = edmParameterSet.getParameter<std::string>("dumpHitsFileName");
      observers.emplace_back(std::make_unique<PatternsPtAssignment>(edmParameterSet, omtfConfig.get(), omtfProcGoldenPat->getPatterns(), ""));
    }*/
  }
}
