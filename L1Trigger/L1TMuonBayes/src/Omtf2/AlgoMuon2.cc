/*
 * AlgoMuon2.cc
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/Omtf2/AlgoMuon2.h"
#include "iomanip"

AlgoMuon2::AlgoMuon2(const ProcConfigurationBase* config): AlgoMuonBase(config) {
  // TODO Auto-generated constructor stub

}

AlgoMuon2::~AlgoMuon2() {
  // TODO Auto-generated destructor stub
}

//TODO fix
void AlgoMuon2::addStubResult(float pdfVal, bool valid, int pdfBin, int layer, MuonStubPtr stub) {
  if(valid) {
    //stubResults.emplace_back(pdfVal, valid, pdfBin, layer, stub);
    displacedLikelihood += pdfVal;
    firedLayerBitsInBx.at(stub->bx).set(layer);
  }
  stubResults[layer] = StubResult(pdfVal, valid, pdfBin, layer, stub);

  //stub result is added even thought it is not valid since this might be needed for debugging or optimization
}

std::ostream & operator<< (std::ostream &out, const AlgoMuon2 &o) {
  out <<"AlgoMuon2: ";
  out << " pt: "   << o.getPtHw()
      << ", phi: " << o.getPhiHw()
      << ", eta: " << o.getEtaHw()
      << ", q: "   << o.getQuality()
//      << ", bx: "  << o.getBx()
      << ", charge: "<< o.getCharge()
      << ", displacement: "  << o.displacementHw
      <<", muonLikelihood "<<o.muonLikelihood
      <<", displacedLikelihood "<<o.displacedLikelihood;

  out <<"\nstubResults: "<<std::endl;
  for(auto& stubResult : o.stubResults) {
    if(stubResult.getMuonStub() ) {
      out <<"layer "<<std::setw(2)<<stubResult.getLayer()<<" valid "<<stubResult.getValid()
          <<" "<<(*stubResult.getMuonStub())<<std::endl;
    }
  }

  return out;
}
