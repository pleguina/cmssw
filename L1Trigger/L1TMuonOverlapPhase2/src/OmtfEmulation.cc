/*
 * OmtfEmulation.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfEmulation.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"

OmtfEmulation::OmtfEmulation(const edm::ParameterSet& edmParameterSet, MuStubsInputTokens& muStubsInputTokens, edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2):
OMTFReconstruction(edmParameterSet, muStubsInputTokens), inputTokenDTPhPhase2(inputTokenDTPhPhase2)
{

}

OmtfEmulation::~OmtfEmulation() {

}

void OmtfEmulation::beginJob() {
  inputMaker.reset(new InputMakerPhase2(edmParameterSet, muStubsInputTokens, inputTokenDTPhPhase2, omtfConfig.get(), new OmtfPhase2AngleConverter() ));
}

void OmtfEmulation::addObservers() {
  OMTFReconstruction::addObservers();

  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>* >(omtfProc.get());
  if(omtfProcGoldenPat) {
/*    if(edmParameterSet.exists("neuralNetworkFile") ) {
      edm::LogImportant("OMTFReconstruction") << "constructing PtAssignmentNN"<<std::endl;
      std::string neuralNetworkFile = edmParameterSet.getParameter<edm::FileInPath>("neuralNetworkFile").fullPath();
      omtfProcGoldenPat->setPtAssignment(new PtAssignmentNN(edmParameterSet, omtfConfig, neuralNetworkFile) ); //TODO change to dynamic_cast and check the type
    }*/

/*    if(edmParameterSet.exists("patternsPtAssignment") && edmParameterSet.getParameter<bool>("patternsPtAssignment")) {
      //std::string rootFileName = edmParameterSet.getParameter<std::string>("dumpHitsFileName");
      observers.emplace_back(std::make_unique<PatternsPtAssignment>(edmParameterSet, omtfConfig.get(), omtfProcGoldenPat->getPatterns(), ""));
    }*/
  }
}
