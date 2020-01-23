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
    std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps, std::string rootFileName):
PatternOptimizerBase(edmCfg, omtfConfig, gps){

initializeTTree(rootFileName);

}

DataROOTDumper2::~DataROOTDumper2() {
  saveTTree();
}


void DataROOTDumper2::initializeTTree(std::string rootFileName) {
  rootFile = new TFile(rootFileName.c_str(), "RECREATE");
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

  ptGenPos = new TH1I("ptGenPos", "ptGenPos", 400, 0, 200);
  ptGenNeg = new TH1I("ptGenNeg", "ptGenNeg", 400, 0, 200);
}

void DataROOTDumper2::saveTTree() {
  rootFile->Write();
  ptGenPos->Write();
  ptGenNeg->Write();

  delete rootFile;
  delete rootTree;
}


void DataROOTDumper2::observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  int muonCharge = 0;
  if(simMuon) {
    if(abs(simMuon->momentum().eta()) < 0.8 || abs(simMuon->momentum().eta()) > 1.24)
      return;

    muonCharge = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;
    if(muonCharge > 0)
      ptGenPos->Fill(simMuon->momentum().pt());
    else
      ptGenNeg->Fill(simMuon->momentum().pt());
  }

  if(simMuon == 0 || !omtfCand->isValid()) //no sim muon or empty candidate
    return;

  //PatternOptimizerBase::observeEventEnd(iEvent, finalCandidates); not needed

  event.muonPt = simMuon->momentum().pt();
  event.muonEta = simMuon->momentum().eta();

/*  if(abs(event.muonEta) < 0.8 || abs(event.muonEta) > 1.24)
    return;*/

  event.muonPhi = simMuon->momentum().phi();
  event.muonCharge = muonCharge;  //TODO


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

        /*int phiDist = stubResult.getPdfBin() - pdfMiddle;

        phiDist = (phiDist << omtfCand->getGoldenPatern()->getDistPhiBitShift(iLogicLayer, omtfCand->getRefLayer()) );
        phiDist += (omtfCand->getGoldenPatern()->meanDistPhiValue(iLogicLayer, omtfCand->getRefLayer()) ); //removing the shift applied in the GoldenPatternBase::process1Layer1RefLayer
        //TODO include the phiDist = phiDist >> this->getDistPhiBitShift(iLayer, iRefLayer); applied in the GoldenPatternBase::process1Layer1RefLaye
        hit.phiDist = phiDist;*/


        int hitPhi = stubResult.getMuonStub()->phiHw;
        unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[omtfCand->getRefLayer()];
        int phiRefHit = gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->phiHw;

        if(omtfConfig->isBendingLayer(iLogicLayer) ) {
          hitPhi = stubResult.getMuonStub()->phiBHw;
          phiRefHit = 0; //phi ref hit for the banding layer set to 0, since it should not be included in the phiDist
        }

        //phiDist = hitPhi - phiRefHit;
        hit.phiDist = hitPhi - phiRefHit;

       /* edm::LogVerbatim("l1tMuBayesEventPrint")<<" muonPt "<<event.muonPt<<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer
            <<" layer "<<int(hit.layer)<<" PdfBin "<<stubResult.getPdfBin()<<" hit.phiDist "<<hit.phiDist<<" valid "<<stubResult.getValid()<<" " //<<" phiDist "<<phiDist
            <<" getDistPhiBitShift "<<omtfCand->getGoldenPatern()->getDistPhiBitShift(iLogicLayer, omtfCand->getRefLayer())
            <<" meanDistPhiValue   "<<omtfCand->getGoldenPatern()->meanDistPhiValue(iLogicLayer, omtfCand->getRefLayer())//<<(phiDist != hit.phiDist? "!!!!!!!<<<<<" : "")
            <<endl;*/

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
