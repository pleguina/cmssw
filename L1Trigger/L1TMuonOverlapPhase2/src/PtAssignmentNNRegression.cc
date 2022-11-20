/*
 * PtAssignmentNN.cc
 *
 *  Created on: May 8, 2020
 *      Author: kbunkow
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <L1Trigger/L1TMuonOverlapPhase2/interface/PtAssignmentNNRegression.h>

#include <sstream>
#include <fstream>

PtAssignmentNNRegression::PtAssignmentNNRegression(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::string networkFile): PtAssignmentBase(omtfConfig) {
  std::ifstream ifs(networkFile);
  edm::LogImportant("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<" networkFile "<<networkFile<<std::endl;
  lutNetworkFP.load(networkFile);
  edm::LogImportant("OMTFReconstruction") <<" "<<__FUNCTION__<<":"<<__LINE__<<std::endl;

  /*
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
*/
}

struct OmtfHit {
  union {
    unsigned long rawData = 0;

    struct {
      char layer;
      char quality;
      char z;
      char valid;
      short eta;
      short phiDist;
    };
  };

  OmtfHit(unsigned long rawData): rawData(rawData) {}
};

bool omtfHitToEventInput(OmtfHit& hit, std::vector<float>& inputs, unsigned int omtfRefLayer, bool print) {
	float offset = (omtfRefLayer<<7) + 64;

	if(hit.valid) {
		if( (hit.layer == 1 || hit.layer == 3 || hit.layer == 5) && hit.quality < 4) ///TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			return false;

		int rangeFactor = 2; //rangeFactor scales the hit.phiDist such that the event->inputs is smaller then 63
		if(hit.layer == 1) {
			rangeFactor = 8;
		}
		/*else if(hit.layer == 8 || hit.layer == 17) {
            rangeFactor = 4;
        }*/
		else if(hit.layer == 3) {
			rangeFactor = 4;
		}
		else if(hit.layer == 9) {
			rangeFactor = 1;
		}
		/*else {
            rangeFactor = 2;
        }
		 */

		rangeFactor *= 2; //TODO !!!!!!!!!!!!!!!!!!!

		if(abs(hit.phiDist) >= (63 * rangeFactor) ) {
			edm::LogImportant("OMTFReconstruction")    //<<" muonPt "<<omtfEvent.muonPt<<" omtfPt "<<omtfEvent.omtfPt
			<<" RefLayer "<<omtfRefLayer<<" layer "
			<<int(hit.layer)<<" hit.phiDist "<<hit.phiDist
			<<" valid "<<((short)hit.valid)<<" !!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			hit.phiDist = copysign(63 * rangeFactor, hit.phiDist);
		}

		inputs.at(hit.layer) = (float)hit.phiDist / (float)rangeFactor + offset;

		if(inputs.at(hit.layer) >= 1022) //the last address i.e. 1023 is reserved for the no-hit value, so interpolation between the 1022 and 1023 has no sense
			inputs.at(hit.layer) = 1022;

		if(print || inputs.at(hit.layer) < 0) {
			edm::LogImportant("OMTFReconstruction") //<<"rawData "<<hex<<setw(16)<<hit.rawData
			<<" layer "<<dec<<int(hit.layer);
			edm::LogImportant("OMTFReconstruction") <<" phiDist "<<hit.phiDist<<" inputVal "<<inputs.at(hit.layer)<<" hit.z "<<int(hit.z)<<" valid "<<((short)hit.valid)<<" quality "<<(short)hit.quality<<" omtfRefLayer "<<omtfRefLayer;
			if(inputs.at(hit.layer) < 0)
				edm::LogImportant("OMTFReconstruction") <<" event->inputs.at(hit.layer) < 0 !!!!!!!!!!!!!!!!!"<<endl;
			edm::LogImportant("OMTFReconstruction") <<endl;
		}

		if(inputs[hit.layer] >= 1024) { //TODO should be the size of the LUT of the first layer
			edm::LogImportant("OMTFReconstruction") <<" event->inputs[hit.layer] >= 1024 !!!!!!!!!!!!!!!!!"<<endl;
		}
		return true;
	}

	return false;

}

PtAssignmentNNRegression::~PtAssignmentNNRegression() {
  // TODO Auto-generated destructor stub
}


std::vector<float> PtAssignmentNNRegression::getPts(const AlgoMuons::value_type& algoMuon) {
  auto& gpResult = algoMuon->getGpResult();
  //int pdfMiddle = 1<<(omtfConfig->nPdfAddrBits()-1);

/*
  edm::LogVerbatim("l1tMuBayesEventPrint")<<"DataROOTDumper2:;observeEventEnd muonPt "<<event.muonPt<<" muonCharge "<<event.muonCharge
      <<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer<<" omtfPtCont "<<event.omtfPtCont
      <<std::endl;
*/

  unsigned int inputCnt = 18; //TDOO!!!!!
  unsigned int outputCnt = 2;
  const float noHitVal = 1023.;

  //edm::LogImportant("OMTFReconstruction") <<"\n----------------------"<<endl;
  //edm::LogImportant("OMTFReconstruction") <<(*algoMuon)<<std::endl;

  std::vector<float> inputs;

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

      OmtfHit hit(0);

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

      omtfHitToEventInput(hit, inputs, algoMuon->getRefLayer(), false);
    }
  }

  std::vector<double> nnResult;
  lutNetworkFP.run(inputs, noHitVal, nnResult);

  double pt = std::copysign(nnResult.at(0), nnResult.at(1));

  std::vector<float> pts;
  pts.emplace_back(pt);
  return pts;

  //event.print();
/*
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

  return pts;*/
}
