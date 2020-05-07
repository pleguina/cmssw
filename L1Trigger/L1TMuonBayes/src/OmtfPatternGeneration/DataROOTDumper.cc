/*
 * DataROOTDumper.cc
 *
 *  Created on: Tue Apr 16 09:57:08 CEST 2019
 *      Author: akalinow
 */
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <boost/range/adaptor/reversed.hpp>
#include <boost/timer/timer.hpp>

#include "TFile.h"
#include "TTree.h"

#include <L1Trigger/L1TMuonBayes/interface/OmtfPatternGeneration/DataROOTDumper.h>

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
DataROOTDumper::DataROOTDumper(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
			       std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  PatternOptimizerBase(edmCfg, omtfConfig, gps){

  initializeTTree();

}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
void DataROOTDumper::initializeTTree(){


  myFile = new TFile("OMTFHits.root","RECREATE");
  myTree = new TTree("OMTFHitsTree","");
  myTree->Branch("muonPt", &myEvent.muonPt);
  myTree->Branch("muonEta", &myEvent.muonEta);
  myTree->Branch("muonPhi", &myEvent.muonPhi);
  myTree->Branch("muonCharge", &myEvent.muonCharge);

  myTree->Branch("omtfPt", &myEvent.omtfPt);
  myTree->Branch("omtfEta", &myEvent.omtfEta);
  myTree->Branch("omtfPhi", &myEvent.omtfPhi);
  myTree->Branch("omtfCharge", &myEvent.omtfCharge);

  myTree->Branch("omtfScore", &myEvent.omtfScore);
  myTree->Branch("omtfQuality", &myEvent.omtfQuality);
  myTree->Branch("omtfRefLayer", &myEvent.omtfRefLayer);
  myTree->Branch("omtfProcessor", &myEvent.omtfProcessor);

  myTree->Branch("hits", &myEvent.hits);
  myTree->Branch("hitsQuality", &myEvent.hitsQuality);
  		
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
void DataROOTDumper::saveTTree(){

  myFile->Write();
  delete myFile;
  delete myTree;
		
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
DataROOTDumper::~DataROOTDumper() {

  saveTTree();

}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
void DataROOTDumper::observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
					      const AlgoMuons& algoCandidates,
					      const AlgoMuons& gbCandidates,
					      const std::vector<l1t::RegionalMuonCand> & candMuons) {
 
  PatternOptimizerBase::observeProcesorEmulation(iProcessor, mtfType, input, algoCandidates, gbCandidates, candMuons);

  std::map<int, int> hwAddrMap = regionalMuonCand.trackAddress();
  std::bitset<21> hitWord(hwAddrMap[0]);
  std::vector<int> tmpPhi, tmpQuality;

  if(omtfCand->getPt()>0 && regionalMuonCand.hwPt()>0 && regionalMuonCand.hwQual()>=12){
    myEvent.omtfPt =  regionalMuonCand.hwPt();
    myEvent.omtfEta = regionalMuonCand.hwEta();
    myEvent.omtfPhi = regionalMuonCand.hwPhi();
    myEvent.omtfCharge = std::pow(-1, regionalMuonCand.hwSign());      
    myEvent.omtfScore = hwAddrMap[2];
    myEvent.omtfQuality = regionalMuonCand.hwQual();
    myEvent.omtfRefLayer  = hwAddrMap[1];
    myEvent.omtfHitsWord = hwAddrMap[0];
    myEvent.omtfProcessor = iProcessor;
    myEvent.hits.clear();
    myEvent.hitsQuality.clear();
    myEvent.pdfResults.clear();
    /*
    myEvent.pdfResults.resize(gPatterns.size());
    int index = 0;
    for(auto& itGP: gPatterns) {
      int ptCode = itGP->key().thePtCode*std::pow(-1,itGP->key().theCharge) + gPatterns.size();     
      unsigned int procIndx = omtfConfig->getProcIndx(iProcessor, mtfType);
      int iRefHit = 0;
      if(!itGP->getResults()[procIndx][iRefHit].isValid()) continue;
      if(itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() < 3) continue;
      int pdfSumTmp = itGP->getResults()[procIndx][iRefHit].getPdfSum();
      myEvent.pdfResults[ptCode] = pdfSumTmp;
    }
    */
    int phi, hwQuality, tmp, sign;
    for(unsigned int iLogicLayer=0;iLogicLayer<omtfConfig->nLayers();++iLogicLayer){
      for(unsigned int iHit=0;iHit<omtfConfig->nInputs();++iHit){
	phi = input.getPhiHw(iLogicLayer, iHit);	
	if(phi!=5400){
	  if(!input.getMuonStub(iLogicLayer, iHit)) hwQuality = input.getMuonStub(iLogicLayer-1, iHit)->qualityHw;
	  else hwQuality = input.getMuonStub(iLogicLayer, iHit)->qualityHw;
	  if(phi<0) {
	    sign = 1;
	    phi = std::abs(phi);
	  }
	  else sign = 0;
	  tmp = (phi<<10) + (iLogicLayer<<5) + iHit;
	  tmp *= std::pow(-1,sign);  
	  myEvent.hits.push_back(tmp);
	  myEvent.hitsQuality.push_back((hwQuality<<10) + (iLogicLayer<<5) + iHit);
	}
      }    
    }
  }
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
void DataROOTDumper::observeEventEnd(const edm::Event& iEvent) {

  //no sim muon or empty candidate
  if(simMuon){
    myEvent.muonPt = simMuon->momentum().pt();
    myEvent.muonEta = simMuon->momentum().eta();
    myEvent.muonPhi = simMuon->momentum().phi();
    myEvent.muonCharge = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;  
  }
  else{
    myEvent.muonPt = -999.0;
    myEvent.muonEta = -999.0;
    myEvent.muonPhi = -999.0;
    myEvent.muonCharge = -999.0;
  }
  myTree->Fill();

  myEvent.omtfPt = 0;
  myEvent.hits.clear();
  myEvent.hitsQuality.clear();
  myEvent.pdfResults.clear();
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////  
void DataROOTDumper::endJob() { }
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////  

