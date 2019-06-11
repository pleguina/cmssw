/*
 * Omtf2inputMaker.cc
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#include <L1Trigger/L1TMuonBayes/interface/Omtf2/Omtf2InputMaker.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

Omtf2InputMaker::Omtf2InputMaker(): OMTFinputMaker() {
  // TODO Auto-generated constructor stub

}

Omtf2InputMaker::~Omtf2InputMaker() {
  // TODO Auto-generated destructor stub
}

void Omtf2InputMaker::addStub(MuonStubPtrs2D& muonStubsInLayers, unsigned int iLayer, unsigned int iInput, MuonStub& stub) {
  unsigned int nMaxMuStubsPerLayer = 8; //TODO add to config
  if(muonStubsInLayers.at(iLayer).size() < nMaxMuStubsPerLayer ) {
    //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" adding stub to layer to iLayer"<<" iLayer "<<stub<<endl;
    muonStubsInLayers.at(iLayer).emplace_back(std::make_shared<const MuonStub>(stub));
    //TODO it looks that sometimes the CSC segments are duplicated, add a handle for this
    /* e.g.
MuonStub:  logicLayer:  7 type:  8 roll: 0 phiHw:   140 phiBHw:    4 etaHw:  -99 etaSigmaHw:   0 qualityHw: 12  bx: 6  timing:  0  CSC   E:2 S:2 R:2 C:9 L:0
MuonStub:  logicLayer:  7 type:  8 roll: 0 phiHw:   140 phiBHw:    4 etaHw:  -99 etaSigmaHw:   0 qualityHw: 12  bx: 6  timing:  0  CSC   E:2 S:2 R:2 C:9 L:0
     */
  }
}
