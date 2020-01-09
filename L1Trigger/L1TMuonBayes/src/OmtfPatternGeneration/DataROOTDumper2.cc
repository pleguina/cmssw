/*
 * DataROOTDumper2.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#include <L1Trigger/L1TMuonBayes/interface/OmtfPatternGeneration/DataROOTDumper2.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <boost/range/adaptor/reversed.hpp>
#include <boost/timer/timer.hpp>
#include "L1Trigger/L1TMuonBayes/interface/OmtfPatternGeneration/DataROOTDumper.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TFile.h"
#include "TTree.h"

DataROOTDumper2::DataROOTDumper2(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
    std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
PatternOptimizerBase(edmCfg, omtfConfig, gps){

initializeTTree();

}

DataROOTDumper2::~DataROOTDumper2() {
  saveTTree();
}


void DataROOTDumper2::initializeTTree() {
  rootFile = new TFile("OMTFHits.root","RECREATE");
  rootTree = new TTree("OMTFHitsTree","");

  rootTree->Branch("muonPt", &event.muonPt);
  rootTree->Branch("muonEta", &event.muonEta);
  rootTree->Branch("muonPhi", &event.muonPhi);
  rootTree->Branch("muonCharge", &event.muonCharge);

  rootTree->Branch("omtfPt", &event.omtfPt);
  rootTree->Branch("omtfEta", &event.omtfEta);
  rootTree->Branch("omtfPhi", &event.omtfPhi);
  rootTree->Branch("omtfCharge", &event.omtfCharge);

  rootTree->Branch("omtfScore", &event.omtfScore);
  rootTree->Branch("omtfQuality", &event.omtfQuality);
  rootTree->Branch("omtfRefLayer", &event.omtfRefLayer);
  rootTree->Branch("omtfProcessor", &event.omtfProcessor);

  rootTree->Branch("hits", &event.hits);
}

void DataROOTDumper2::saveTTree() {
  rootFile->Write();
  delete rootFile;
  delete rootTree;
}


void DataROOTDumper2::observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  if(simMuon == 0 || !omtfCand->isValid()) //no sim muon or empty candidate
    return;

  //PatternOptimizerBase::observeEventEnd(iEvent, finalCandidates); not needed

  event.muonPt = simMuon->momentum().pt();
  event.muonEta = simMuon->momentum().eta();

  if(abs(event.muonEta) < 0.8 || abs(event.muonEta) > 1.24)
    return;

  event.muonPhi = simMuon->momentum().phi();
  event.muonCharge = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;  //TODO


  if(omtfCand->getPt() > 0) { //&& omtfCand->getFiredLayerCnt() > 3
    event.omtfPt = (omtfCand->getPt() -1)/2.0; //TODO check
    event.omtfEta = omtfCand->getEtaHw()/240.*2.61;
    event.omtfPhi = omtfCand->getPhi();
    event.omtfCharge = std::pow(-1, regionalMuonCand.hwSign());
    event.omtfScore = omtfCand->getPdfSum();
    event.omtfQuality = omtfCand->getQ();
    event.omtfRefLayer  = omtfCand->getRefLayer();
    event.omtfHitsWord = omtfCand->getFiredLayerBits();
    event.omtfProcessor = candProcIndx;

    event.hits.clear();

    //unsigned int iRefHit = omtfCand->getRefHitNumber();

    auto& gpResult = omtfCand->getGpResult();
    int pdfMiddle = 1<<(omtfConfig->nPdfAddrBits()-1);

    for(unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
      auto& stubResult = gpResult.getStubResults()[iLogicLayer];
      if(stubResult.getMuonStub() && stubResult.getValid()) {
        OmtfEvent::Hit hit;
        hit.layer = iLogicLayer;
        hit.quality = stubResult.getMuonStub()->qualityHw;
        hit.eta = stubResult.getMuonStub()->etaHw; //in which scale?
        //hit.phi = gpResult.stubResults[iLogicLayer].
          //  int phiDist
        int phiDist = gpResult.getStubResults()[iLogicLayer].getPdfBin();
        phiDist += omtfCand->getGoldenPatern()->meanDistPhiValue(iLogicLayer, omtfCand->getRefLayer()) - pdfMiddle; //removing the shift applied in the GoldenPatternBase::process1Layer1RefLayer
        //TODO include the phiDist = phiDist >> this->getDistPhiBitShift(iLayer, iRefLayer); applied in the GoldenPatternBase::process1Layer1RefLaye
        hit.phiDist = phiDist;

        if(hit.phiDist > 504 || hit.phiDist < -512 ) {
          edm::LogVerbatim("l1tMuBayesEventPrint")<<" muonPt "<<event.muonPt<<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer
              <<" layer "<<int(hit.layer)<<" hit.phiDist "<<hit.phiDist<<" valid "<<stubResult.getValid()<<" !!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

        }

        stubResult.getMuonStub()->detId;

        DetId detId(stubResult.getMuonStub()->detId);
        if(detId.subdetId() == MuonSubdetId::CSC) {
          CSCDetId cscId(detId);
          hit.z = cscId.chamber() % 2;
        }

        event.hits.push_back(hit.rawData);
      }
    }
    rootTree->Fill();
    evntCnt++;
  }
}

void DataROOTDumper2::endJob() {
  edm::LogVerbatim("l1tMuBayesEventPrint")<<" evntCnt "<<evntCnt<<endl;
}
