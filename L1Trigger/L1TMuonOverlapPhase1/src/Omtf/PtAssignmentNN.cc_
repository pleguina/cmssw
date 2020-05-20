/*
 * PtAssignmentNN.cc
 *
 *  Created on: May 8, 2020
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/Omtf/PtAssignmentNN.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "lutNN/lutNN2/interface/NetworkSerialization.h"
#include "lutNN/lutNN2/interface/EventsGeneratorOmtf.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <sstream>
#include <fstream>

PtAssignmentNN::PtAssignmentNN(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::string networkFile): PtAssignmentBase(omtfConfig), network(nullptr) {
  std::ifstream ifs(networkFile);
  {
    edm::LogImportant("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<" networkFile "<<networkFile<<std::endl;
      boost::archive::text_iarchive ia(ifs);
      lutNN::registerClasses(ia);
      ptBins.clear();
      ia >> ptBins;
      edm::LogImportant("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<std::endl;
      ia >> network;
      edm::LogImportant("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<std::endl;
      // archive and stream closed when destructors are called
  }


  if(edmCfg.exists("nn_pThresholds") ) {
    auto nn_pThresholds = edmCfg.getParameter<vector<double> >("nn_pThresholds");

    TFile* ptCalibrationFile =  nullptr;
    if(edmCfg.exists("ptCalibrationFileName") ) {
      auto ptCalibrationFileName = edmCfg.getParameter<edm::FileInPath>("ptCalibrationFileName").fullPath();
      ptCalibrationFile =  new TFile(ptCalibrationFileName.c_str());
      edm::LogImportant("OMTFReconstruction") <<" using ptCalibrationFile "<<ptCalibrationFileName<<std::endl;
    }

    if(ptCalibrationFile ==  nullptr) {
      edm::LogImportant("OMTFReconstruction") <<" no ptCalibrationFile!!!!!!"<<std::endl;
    }

    classifierToRegressions.reserve(nn_pThresholds.size() );
    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream name;
      name<<"ptCalibrationPSumInter"<<nn_pThreshold;
      std::string nameStr = boost::replace_all_copy(name.str(), ".", "");
      lutNN::PtCalibration* ptCalibration = new lutNN::PtCalibration(ptCalibrationFile, nameStr);

      classifierToRegressions.emplace_back(make_unique<lutNN::ClassifierToRegressionPSum>(ptBins, nn_pThreshold, ptCalibration));
      edm::LogImportant("OMTFReconstruction") <<" adding classifierToRegression, nn_pThreshold "<<nn_pThreshold<<std::endl;
    }
  }
  edm::LogImportant("OMTFReconstruction") <<" PtAssignmentNN created, networkFile "<<networkFile<<" classifierToRegressions.size() "<<classifierToRegressions.size()<<std::endl;
}

PtAssignmentNN::~PtAssignmentNN() {
  // TODO Auto-generated destructor stub
}


std::vector<float> PtAssignmentNN::getPts(const AlgoMuons::value_type& algoMuon) {
  auto& gpResult = algoMuon->getGpResult();
  //int pdfMiddle = 1<<(omtfConfig->nPdfAddrBits()-1);

/*
  edm::LogVerbatim("l1tMuBayesEventPrint")<<"DataROOTDumper2:;observeEventEnd muonPt "<<event.muonPt<<" muonCharge "<<event.muonCharge
      <<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer<<" omtfPtCont "<<event.omtfPtCont
      <<std::endl;
*/

  unsigned int inputCnt = 18; //TDOO!!!!!
  unsigned int outputCnt = network.getOutputNode()->getOutputValues().size();
  const float noHitVal = 1023.;

  lutNN::EventFloat event(inputCnt, outputCnt, noHitVal);

  //cout<<"\n----------------------"<<endl;
  //cout<<(*algoMuon)<<std::endl;

  for(unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
    auto& stubResult = gpResult.getStubResults()[iLogicLayer];
    if(stubResult.getMuonStub() ) {//&& stubResult.getValid() //TODO!!!!!!!!!!!!!!!!1
      int hitPhi = stubResult.getMuonStub()->phiHw;
      unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[algoMuon->getRefLayer()];
      int phiRefHit = gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->phiHw;

      if(omtfConfig->isBendingLayer(iLogicLayer) ) {
        hitPhi = stubResult.getMuonStub()->phiBHw;
        phiRefHit = 0; //phi ref hit for the banding layer set to 0, since it should not be included in the phiDist
      }

      lutNN::OmtfEvent::Hit hit(0);

      hit.layer = iLogicLayer;
      hit.quality = stubResult.getMuonStub()->qualityHw;
      hit.eta = stubResult.getMuonStub()->etaHw; //in which scale?
      hit.valid = stubResult.getValid();

      //phiDist = hitPhi - phiRefHit;
      hit.phiDist = hitPhi - phiRefHit;

     /* edm::LogVerbatim("l1tMuBayesEventPrint")<<" muonPt "<<event.muonPt<<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer
          <<" layer "<<int(hit.layer)<<" PdfBin "<<stubResult.getPdfBin()<<" hit.phiDist "<<hit.phiDist<<" valid "<<stubResult.getValid()<<" " //<<" phiDist "<<phiDist
          <<" getDistPhiBitShift "<<omtfCand->getGoldenPatern()->getDistPhiBitShift(iLogicLayer, omtfCand->getRefLayer())
          <<" meanDistPhiValue   "<<omtfCand->getGoldenPatern()->meanDistPhiValue(iLogicLayer, omtfCand->getRefLayer())//<<(phiDist != hit.phiDist? "!!!!!!!<<<<<" : "")
          <<endl;*/

      if(hit.phiDist > 504 || hit.phiDist < -512 ) {
        edm::LogVerbatim("l1tMuBayesEventPrint")
        //<<" muonPt "<<event.muonPt<<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer
            <<" layer "<<int(hit.layer)<<" hit.phiDist "<<hit.phiDist<<" valid "<<stubResult.getValid()<<" !!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      }

      stubResult.getMuonStub()->detId;

      DetId detId(stubResult.getMuonStub()->detId);
      if(detId.subdetId() == MuonSubdetId::CSC) {
        CSCDetId cscId(detId);
        hit.z = cscId.chamber() % 2;
      }

      lutNN::EventsGeneratorOmtf::omtfHitToEventInput(hit, &event, algoMuon->getRefLayer(), false);
    }
  }

  network.run(&event);

  //event.print();

  std::vector<float> pts(classifierToRegressions.size(), 0);

  unsigned int i =0;
  for(auto& classifierToRegression : classifierToRegressions) {
    auto orgValue = classifierToRegression->getValue(&event);
    auto absOrgValue = fabs(orgValue);
    pts.at(i) = classifierToRegression->getCalibratedValue(absOrgValue);
    pts.at(i) = std::copysign(pts.at(i), orgValue);

    LogTrace("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<" orgValue "<<orgValue<<" pts["<<i<<"] "<<pts[i]<<std::endl;
    //std::cout<<"nn pts["<<i<<"] "<<pts[i]<< std::endl;
    i++;
  }

  return pts;
}
