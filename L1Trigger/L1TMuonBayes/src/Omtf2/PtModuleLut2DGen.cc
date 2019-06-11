/*
 * PtModuleLut2DGen.cc
 *
 *  Created on: May 22, 2019
 *      Author: kbunkow
 */

#include <L1Trigger/L1TMuonBayes/interface/Omtf2/PtModuleLut2DGen.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

PtModuleLut2DGen::PtModuleLut2DGen(const ProcConfigurationBase* config): PtModuleLut2D(config) {

}

PtModuleLut2DGen::~PtModuleLut2DGen() {
  // TODO Auto-generated destructor stub
}

void PtModuleLut2DGen::prepare(Mode mode, SampleType sampleType) {
  //adding one one more LUT val to  keep number of events
  this->mode =  mode;
  this->sampleType = sampleType;

  for(auto& lut : lut2ds) {
    if(lut->getLutValues().size() < 5) {
      for(auto& lut : lut2ds) {
        lut->getLutValues().emplace_back(std::vector<std::vector<Lut2dInter::lutValueType> >( 1<<lut->getLutConfig().addr1.addrBits, std::vector<Lut2dInter::lutValueType>(1<<lut->getLutConfig().addr2.addrBits, 0) ));
      }
    }

    for(unsigned int iOut = 0; iOut < lut->getLutValues().size(); iOut++) {
      for(unsigned int addr1 = 0; addr1 < lut->getLutValues()[iOut].size(); addr1++) {
        for(unsigned int addr2 = 0; addr2 < lut->getLutValues()[iOut][addr1].size(); addr2++) {
          if(mode == Mode::calcualteMeans) {
            if(iOut == 0 || iOut == 2)
              lut->getLutValues()[iOut][addr1][addr2] = 0;
          }
          else if(mode == Mode::calcualteStdDevs) {
            if(iOut == 1 || iOut == 3)
              lut->getLutValues()[iOut][addr1][addr2] = 0;
          }

          if(iOut == 4 &&  sampleType == SampleType::muons)
            lut->getLutValues()[iOut][addr1][addr2] = 0;


          if(iOut == 5 &&  sampleType == SampleType::displacedMuons)
            lut->getLutValues()[iOut][addr1][addr2] = 0;
        }
      }
    }
  }
}

AlgoMuon2Ptr PtModuleLut2DGen::processStubs(const MuonStubsInput& muonStubs) {
  AlgoMuon2Ptr algoMuon = std::make_shared<AlgoMuon2>(config);
/*
  if(!simMuon)
    return algoMuon;

  int charge = 0;
  if(simMuon->type() == 13) //muon
    charge = -1;
  else if(simMuon->type() == -13) //muon
    charge = 1;

  double muPt = simMuon->momentum().Pt() * charge; //TODO convert pt
*/
  if(!genMuon)
    return algoMuon;

  int charge = 0;
  if(genMuon->pdgId() == 13) //muon
    charge = -1;
  else if(genMuon->pdgId() == -13) //muon
    charge = 1;

  double muPt = charge / genMuon->pt(); //TODO convert pt

  double muDxy = (-1 * genMuon->vx() * genMuon->py() + genMuon->vy() * genMuon->px()) / genMuon->pt();;

  //double muVt = sqrt( genMuon->vx() * genMuon->vx() + genMuon->vy() * genMuon->vy()); //vertex distance from the beam

  double eventWeigt = 1;
  //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" muPt "<<muPt<<" "<<std::endl;
  for(auto& lut : lut2ds) {
    if(mode == Mode::calcualteMeans) {
      lut->updateValue(muonStubs, 0, muPt * eventWeigt);

      muDxy = round(muDxy); //TDOO remove
      lut->updateValue(muonStubs, 2, muDxy * eventWeigt);
/*
      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" "<<(*lut)<<" === result: pt "<<std::setw(10)<<outVals.at(0)<<" pt weight "<<std::setw(10)<<outVals.at(1)
                                                              <<" || displacement "<<std::setw(10)<<outVals.at(2)<<" displWeight "<<std::setw(10)<<outVals.at(3)<<std::endl;
*/
    }
    else if(mode == Mode::calcualteStdDevs) {
      std::vector<float> outVals = lut->getValues(muonStubs);
      assert(outVals.size() == 6);

      double ptStdDev = pow(outVals[0] - muPt, 2);
      double dispStdDev = pow(outVals[2] - muDxy, 2);

      bool upadatedPt = lut->updateValue(muonStubs, 1, ptStdDev * eventWeigt);

      lut->updateValue(muonStubs, 3, dispStdDev * eventWeigt);

      if(upadatedPt)
        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" "<<(*lut)<<" === result: pt "<<std::setw(10)<<outVals.at(0)<<" ptStdDev "<<std::setw(10)<<ptStdDev<<" muPt "<<std::setw(10)<<muPt
                                                              //<<" || displacement "<<std::setw(10)<<outVals.at(2)<<" displWeight "<<std::setw(10)<<outVals.at(3)
                                                              <<std::endl;

    }

    if(sampleType == SampleType::muons)
      lut->updateValue(muonStubs, 4, eventWeigt);


    if(sampleType == SampleType::displacedMuons)
      lut->updateValue(muonStubs, 5, eventWeigt);

    //LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  }

  return algoMuon;
}

void PtModuleLut2DGen::endJob() {
  edm::LogVerbatim("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<std::endl;

  int minEventCnt = 1; // will set values to 0 when calculating logLikelihood

  for(auto& lut : lut2ds) {
    for(unsigned int iOut = 0; iOut < lut->getLutValues().size() -1; iOut++) {
      if(mode == Mode::calcualteMeans) {
        if(iOut != 0 && iOut != 2)
          continue;
      }
      else if(mode == Mode::calcualteStdDevs) {
        if(iOut != 1 && iOut != 3)
          continue;
      }

      for(unsigned int addr1 = 0; addr1 < lut->getLutValues()[iOut].size(); addr1++) {
        for(unsigned int addr2 = 0; addr2 < lut->getLutValues()[iOut][addr1].size(); addr2++) {
          double eventCnt = lut->getLutValues().at(4)[addr1][addr2];

          if(iOut == 0 || iOut == 2) {//mean pt or disp
            double val = lut->getLutValues()[iOut][addr1][addr2];
            if(eventCnt >= minEventCnt) {
              lut->getLutValues()[iOut][addr1][addr2] = val /eventCnt;
            }
            else {
              lut->getLutValues()[iOut][addr1][addr2] = 0; //TODO is it ok?
            }

            edm::LogVerbatim("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" calcualteMeans "<<(*lut)
                          <<" iOut "<<" addr1 "<<addr1<<" addr2 "<<addr2<<" val "<<val<<" eventCnt "<<eventCnt<<" lutVal "<<lut->getLutValues()[iOut][addr1][addr2]<<std::endl;
          }
          else if(iOut == 1 || iOut == 3) {//weights i.e. 1/stdDev
            double val = lut->getLutValues()[iOut][addr1][addr2];
            if(eventCnt < minEventCnt) {
              lut->getLutValues()[iOut][addr1][addr2] = 0;
            }
            else if(val == 0) { //StdDev is zero!!
              lut->getLutValues()[iOut][addr1][addr2] = 10000;//0xffff; //TODO put some reasonable value
            }
            else {
              lut->getLutValues()[iOut][addr1][addr2] = eventCnt / val;
            }

            edm::LogVerbatim("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" calcualteStdDevs"<<(*lut)
                                      <<" iOut "<<" addr1 "<<addr1<<" addr2 "<<addr2<<" val "<<val<<" eventCnt "<<eventCnt<<" lutVal "<<lut->getLutValues()[iOut][addr1][addr2]<<std::endl;

          }
        }
      }
    }

  }

  if(mode == Mode::calcualteStdDevs) {
    if(sampleType == SampleType::muons) {
      calcualteLogLikelihoods(4);
    }
    else if(sampleType == SampleType::displacedMuons) {
      calcualteLogLikelihoods(5);
    }
  }
}


void PtModuleLut2DGen::calcualteLogLikelihoods(unsigned int iOut) {
  unsigned int  pdfMaxLogVal = 511; //TODO optimize
  double minPdfVal = 0.00003;      //TODO optimize

  for(auto& lut : lut2ds) {
    double norm = 0;
      for(unsigned int addr1 = 0; addr1 < lut->getLutValues().at(iOut).size(); addr1++) {
        for(unsigned int addr2 = 0; addr2 < lut->getLutValues().at(iOut)[addr1].size(); addr2++) {
          norm += lut->getLutValues().at(iOut)[addr1][addr2];
        }
      }
      for(unsigned int addr1 = 0; addr1 < lut->getLutValues().at(iOut).size(); addr1++) {
        for(unsigned int addr2 = 0; addr2 < lut->getLutValues().at(iOut)[addr1].size(); addr2++) {
          int eventCnt = lut->getLutValues().at(iOut)[addr1][addr2];
          double pdfVal = lut->getLutValues().at(iOut)[addr1][addr2] / norm;

          const double minPlog =  log(minPdfVal);
          //const double pdfMaxLogVal = config->pdfMaxLogValue();

          double logPdf = 0;
          if(pdfVal >= minPdfVal)  //for removing points with small statistics - TODO tune
            //if(pdfVal > 0)
          {
            logPdf = pdfMaxLogVal - log(pdfVal) / minPlog * pdfMaxLogVal;

            lut->getLutValues().at(iOut)[addr1][addr2] = round(logPdf); //round(100. * lut->getLutValues().at(iOut)[addr1][addr2] / norm); //TODO change to log likelihood

            edm::LogVerbatim("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" iOut "<<iOut <<(*lut)
                                                 <<" iOut "<<" addr1 "<<addr1<<" addr2 "<<addr2<<" eventCnt "<<eventCnt<<" logPdf "<<logPdf<<" lutValue "<<lut->getLutValues()[iOut][addr1][addr2]<<std::endl;
          }
          else {
            for(unsigned int iOut1 = 0; iOut1 < lut->getLutValues().size(); iOut1++) {
              lut->getLutValues().at(iOut1)[addr1][addr2] = 0;
            }

            edm::LogVerbatim("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" iOut "<<iOut <<(*lut)
                                                             <<" iOut "<<" addr1 "<<addr1<<" addr2 "<<addr2<<" eventCnt "<<eventCnt<<" logPdf "<<logPdf<<" lutValue "<<lut->getLutValues().at(iOut)[addr1][addr2]<<std::endl;
          }
        }
      }
  }
}
