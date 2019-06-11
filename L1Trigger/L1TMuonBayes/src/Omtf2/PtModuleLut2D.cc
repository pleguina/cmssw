/*
 * TrackMeasureModule.cc
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/Omtf2/PtModuleLut2D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <memory>
#include <iomanip>

#include <TH2F.h>

#undef BOOST_DISABLE_ASSERTS //TODO remove for production version

PtModuleLut2D::PtModuleLut2D(const ProcConfigurationBase* config): config(config) {
  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" config "<<config<<std::endl;
}

PtModuleLut2D::~PtModuleLut2D() {
  // TODO Auto-generated destructor stub
}

void PtModuleLut2D::init() {
  Lut2dInter::LutConfig lc;

  lc.outValCnt = 6; //the value 4 is for storing the number of events (lieklihhod) for normal muon, the value 5 is for number of events displaced muons
  lc.addr1.addrBits = 5;
  lc.addr1.restBits = 5;
  lc.addr2.addrBits = 5;
  lc.addr2.restBits = 5;


  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 0;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 2;
  lc.addr2.bitShift = 2;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 4;
  lc.addr2.bitShift = 2;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 3;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 2;
  lc.addr2.layerB = 4;
  lc.addr2.bitShift = 3;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 0;

  lc.addr2.layerA = 3;
  lc.addr2.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr2.bitShift = 1;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 5;
  lc.addr2.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr2.bitShift = 1;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 0;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 6;
  lc.addr2.bitShift = 2;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 7;
  lc.addr2.bitShift = 3;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 8;
  lc.addr2.bitShift = 3;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 1;
  lc.addr1.layerB = Lut2d::LutConfig::NOT_A_LAYER;
  lc.addr1.bitShift = 1;

  lc.addr2.layerA = 0;
  lc.addr2.layerB = 14;
  lc.addr2.bitShift = 2;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 0;
  lc.addr1.layerB = 2;
  lc.addr1.bitShift = 3;

  lc.addr2.layerA = 2;
  lc.addr2.layerB = 4;
  lc.addr2.bitShift = 3;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 0;
  lc.addr1.layerB = 2;
  lc.addr1.bitShift = 3;

  lc.addr2.layerA = 2;
  lc.addr2.layerB = 6;
  lc.addr2.bitShift = 3;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));



  lc.addr1.layerA = 6;
  lc.addr1.layerB = 7;
  lc.addr1.bitShift = 5;

  lc.addr2.layerA = 7;
  lc.addr2.layerB = 8;
  lc.addr2.bitShift = 4;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 9;
  lc.addr1.layerB = 7;
  lc.addr1.bitShift = 5;

  lc.addr2.layerA = 7;
  lc.addr2.layerB = 8;
  lc.addr2.bitShift = 4;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 16;
  lc.addr1.layerB = 7;
  lc.addr1.bitShift = 5;

  lc.addr2.layerA = 7;
  lc.addr2.layerB = 8;
  lc.addr2.bitShift = 4;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 7;
  lc.addr1.layerB = 8;
  lc.addr1.bitShift = 5;

  lc.addr2.layerA = 8;
  lc.addr2.layerB = 17;
  lc.addr2.bitShift = 5;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 10;
  lc.addr1.layerB = 11;
  lc.addr1.bitShift = 4;

  lc.addr2.layerA = 11;
  lc.addr2.layerB = 2;
  lc.addr2.bitShift = 2;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));


  lc.addr1.layerA = 10;
  lc.addr1.layerB = 2;
  lc.addr1.bitShift = 2;

  lc.addr2.layerA = 3;
  lc.addr2.layerB = Lut2d::LutConfig::NOT_A_LAYER;;
  lc.addr2.bitShift = 1;
  lut2ds.emplace_back(std::make_shared<Lut2dInter>(lc));

}

//each bit in the interpolate switch on (if 1) or off (if 0) the interpolation in calculation each output value from the LUT
void PtModuleLut2D::setInterpolate(const boost::dynamic_bitset<>& interpolate) {
  for(auto& lut : lut2ds) {
    lut->setInterpolate(interpolate);
  }
}

AlgoMuon2Ptr PtModuleLut2D::processStubs(const MuonStubsInput& muonStubs) {
  AlgoMuon2Ptr algoMuon = std::make_shared<AlgoMuon2>(config);

  float pt = 0;
  float ptWeightSum = 0;

  float displacement = 0;
  float displacementWeightSum = 0;

  float muonLikelihood = 0;
  float displLikelihood = 0;

  for(auto& lut : lut2ds) {
    std::vector<float> outVals = lut->getValues(muonStubs);
    assert(outVals.size() >= 4);

    if(weightPt) {
      pt += outVals[0] * outVals[1]; //outVals[0] = pt, outVals[1] = weight
      ptWeightSum += outVals[1];
    }
    else { //dont use weights
      pt += outVals[0]; //outVals[0] = pt, outVals[1] = weight
      ptWeightSum += 1;
    }

    if(weightDisplacement) {
      displacement += outVals[2] * outVals[3]; //outVals[0] = pt, outVals[1] = weight
      displacementWeightSum += outVals[3];
    }
    else {
      displacement += outVals[2]; //outVals[0] = pt, outVals[1] = weight
      displacementWeightSum += 1;
    }

    muonLikelihood += outVals[4];
    displLikelihood += outVals[5];

    if(outVals[1] != 0) //if ptWeight != 0 then there were some good stubs
      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" "<<(*lut)<<" === result: pt "<<std::setw(10)<<1./outVals.at(0)<<" pt weight "<<std::setw(10)<<outVals.at(1)
                                                          <<" || displacement "<<std::setw(10)<<outVals.at(2)<<" displWeight "<<std::setw(10)<<outVals.at(3)
                                                          <<" muonLikelihood "<<std::setw(10)<<outVals.at(4)<<" displLikelihood "<<std::setw(10)<<outVals.at(5)<<std::endl;
  }

  if(ptWeightSum)
  {
    //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" pt "<<pt<<" ptWeightSum "<<ptWeightSum<<" pt / ptWeightSum "<<pt / ptWeightSum<<std::endl;
    pt = pt / ptWeightSum;
    //pt = round(pt); changes to always 0 because it pt is less ten 1, TODO multiply by something
  }

  if(displacementWeightSum) {
    displacement = displacement / displacementWeightSum; //todo do it like in firmware
    displacement = round(displacement);
  }

  if(ptWeightSum || displacementWeightSum) {//if ptWeightSum = 0 then it means then none of the LUTs sound anything reasonable, TODO use likelihood instead
    algoMuon->setPtHw( 1./pt); //config->ptGevToHw(abs() ) TODO

    algoMuon->setCharge( pt > 0 ? 1 : -1);

    algoMuon->setDisplacementHw(displacement);
    algoMuon->setValid(true);

    algoMuon->setMuonLikelihood(muonLikelihood);
    algoMuon->setDisplacedLikelihood(displLikelihood);

    LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" pt "<<1./pt<<" ptWeightSum "<<ptWeightSum<<std::endl;
    LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" found "<<(*algoMuon)<<std::endl;
  }

  //set phi, eta,

  return algoMuon;
}

void PtModuleLut2D::saveAsRootHists() {
  edm::Service<TFileService> fs;

  std::vector<std::string> layerNames = {
      "MB1"    ,  //0
      "MB1b"   ,  //1
      "MB2"    ,  //2
      "MB2b"   ,  //3
      "MB3"    ,  //4
      "MB3b"   ,  //5
      "ME1/3"  ,  //6
      "ME2/2"  ,  //7
      "ME3/2"  ,  //8
      "ME1/2"  ,  //9
      "RB1in"  ,  //10
      "RB1out" ,  //11
      "RB2in"  ,  //12
      "RB2out" ,  //13
      "RB3"    ,  //14
      "RE1/3"  ,  //15
      "RE2/3"  ,  //16
      "RE3/3"     //17
  };


  TFileDirectory subDir = fs->mkdir("PtModuleLut2D");
  for(auto& lut : lut2ds) {
    for(unsigned int iOut = 0; iOut < lut->getLutValues().size(); iOut++) {
      std::ostringstream ostr;
      if(iOut == 0)
        ostr<<lut->getName()<<"_pt";
      else if(iOut == 1)
        ostr<<lut->getName()<<"_ptWeight";
      else if(iOut == 2)
        ostr<<lut->getName()<<"_disp";
      else if(iOut == 3)
        ostr<<lut->getName()<<"_disptWeight";
      else if(iOut == 4)
        ostr<<lut->getName()<<"_eventCnt_muons";
      else if(iOut == 5)
        ostr<<lut->getName()<<"_eventCnt_dispMuons";

      TH2F* hist = subDir.make<TH2F>(ostr.str().c_str(), ostr.str().c_str(),
          lut->getLutConfig().addr1.getBinsCnt(), -.5, lut->getLutConfig().addr1.getBinsCnt()-0.5,
          lut->getLutConfig().addr2.getBinsCnt(), -.5, lut->getLutConfig().addr2.getBinsCnt()-0.5 );

      std::string axisName = layerNames[lut->getLutConfig().addr1.layerA];
      if(lut->getLutConfig().addr1.layerB != Lut2d::LutConfig::NOT_A_LAYER) {
        axisName = axisName + " - " + layerNames[lut->getLutConfig().addr1.layerB];
      }
      hist->GetXaxis()->SetTitle(axisName.c_str());

      axisName = layerNames[lut->getLutConfig().addr2.layerA];
      if(lut->getLutConfig().addr2.layerB != Lut2d::LutConfig::NOT_A_LAYER) {
        axisName = axisName + " - " + layerNames[lut->getLutConfig().addr2.layerB];
      }
      hist->GetYaxis()->SetTitle(axisName.c_str());

      for(unsigned int addr1 = 0; addr1 < lut->getLutValues()[iOut].size(); addr1++) {
        for(unsigned int addr2 = 0; addr2 < lut->getLutValues()[iOut][addr1].size(); addr2++) {
          hist->Fill(addr1, addr2, lut->getLutValues()[iOut][addr1][addr2]);
        }
      }
    }
  }
}
