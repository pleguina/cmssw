/*
 * PatternOptimizer.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/PatternOptimizer.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"

#include "Math/VectorUtil.h"

PatternOptimizer::PatternOptimizer(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  edmCfg(edmCfg), omtfConfig(omtfConfig), goldenPatterns(gps), simMuon(0),
  //TODO set desire function here, see https://www.cprogramming.com/c++11/c++11-lambda-closures.html
  updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatVoter_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtDiff_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtDiff2_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtLogDiff_2(omtfCandGp, exptCandGp); } ),
  //updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_1(gp, iLayer, iRefLayer, learingRate); })
  //updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_2(gp, iLayer, iRefLayer, learingRate); })
  updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsVoter_1(gp, iLayer, iRefLayer, learingRate); })
{

}

PatternOptimizer::~PatternOptimizer() {
}

void PatternOptimizer::observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
    const std::vector<AlgoMuon>& algoCandidates,
    std::vector<AlgoMuon>& gbCandidates,
    const std::vector<l1t::RegionalMuonCand> & candMuons) {

  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;
  int patNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(patNum).get(); // expected pattern

  //bool found = false;
  for(auto& gbCandidate : gbCandidates) {
    int refHitNum = gbCandidate.getRefHitNumber();
    if(gbCandidate.getGoldenPatern() != 0 && gbCandidate.getGoldenPatern()->getResults().at(refHitNum).getFiredLayerCnt() > omtfResult.getFiredLayerCnt() ) {
      //cout<<__FUNCTION__<<":"<<__LINE__<<" gbCandidate "<<gbCandidate<<" "<<std::endl;
      omtfCand = gbCandidate;
      omtfResult = gbCandidate.getGoldenPatern()->getResults().at(refHitNum);
      exptResult = exptCandGp->getResults().at(refHitNum);
      //found = true;
    }
  }

  //////////////////////debug printout/////////////////////////////
  /*if(found) {
    GoldenPatternWithStat* omtfCandGp = static_cast<GoldenPatternWithStat*>(omtfCand.getGoldenPatern());
    if( omtfCandGp->key().thePt > 100 && exptCandGp->key().thePt <= 15 ) {
      //std::cout<<iEvent.id()<<std::endl;
      cout<<" ptSim "<<ptSim<<" chargeSim "<<chargeSim<<std::endl;
      std::cout<<"iProcessor "<<iProcessor<<" exptCandGp "<<exptCandGp->key()<<std::endl;
      std::cout<<"iProcessor "<<iProcessor<<" omtfCandGp "<<omtfCandGp->key()<<std::endl;
      std::cout<<"omtfResult "<<std::endl<<omtfResult<<std::endl;
      int refHitNum = omtfCand.getRefHitNumber();
      std::cout<<"other gps results"<<endl;
      for(auto& gp : goldenPatterns) {
        if(omtfResult.getFiredLayerCnt() == gp->getResults().at(refHitNum).getFiredLayerCnt() )
        {
          cout<<gp->key()<<std::endl<<gp->getResults().at(refHitNum)<<std::endl;
        }
      }
      std::cout<<std::endl;
    }
  }*/
}

void PatternOptimizer::observeEventBegin(const edm::Event& iEvent) {
  omtfCand = AlgoMuon();
  omtfResult =  GoldenPatternResult();
  exptResult =  GoldenPatternResult();

  simMuon = findSimMuon(iEvent);
  //cout<<__FUNCTION__<<":"<<__LINE__<<" simMuon "<<simMuon->momentum().pt()<<" "<<std::endl;
}

void PatternOptimizer::observeEventEnd(const edm::Event& iEvent) {
  if(simMuon == 0 || omtfCand.getGoldenPatern() == 0)//no sim muon or empty candidate
    return;

  //TODO move the eta cut somewhere else
/*  if( abs(simMuon->momentum().eta()) < 0.85 || abs(simMuon->momentum().eta()) > 1.25 )
    return;*/

  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCand "<<omtfCand<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfResult "<<std::endl<<omtfResult<<std::endl;

  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;

  GoldenPatternWithStat* omtfCandGp = static_cast<GoldenPatternWithStat*>(omtfCand.getGoldenPatern());

  int patNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(patNum).get(); // expected pattern

  updateStatFunc(omtfCandGp, exptCandGp);

  if( omtfCandGp->key().thePt > 40 && exptCandGp->key().thePt <= 15 ) {
    std::cout<<iEvent.id()<<std::endl;
    std::cout<<" ptSim "<<ptSim<<" chargeSim "<<chargeSim<<" patNum "<<patNum<<std::endl;
    std::cout<<"exptCandGp "<<exptCandGp->key()<<std::endl;
    std::cout<<"exptResult "<<std::endl<<exptResult<<std::endl;

    std::cout<<"omtfCandGp "<<omtfCandGp->key()<<std::endl;
    std::cout<<"omtfResult "<<std::endl<<omtfResult<<std::endl;
    int refHitNum = omtfCand.getRefHitNumber();
    /*std::cout<<"other gps results"<<endl; will not work here , since the gps hold now results for other processor
    for(auto& gp : goldenPatterns) {
      if(omtfResult.getFiredLayerCnt() == gp->getResults().at(refHitNum).getFiredLayerCnt() )
      {
        cout<<gp->key()<<std::endl<<gp->getResults().at(refHitNum)<<std::endl;
      }
    }*/
    std::cout<<std::endl;
  }

}

void PatternOptimizer::endJob() {
  double learingRate = 0.1;
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        updatePdfsFunc(gp.get(), iLayer, iRefLayer, learingRate);
      }
    }
  }

  //"optimisedPats.xml";
  std::string fName = edmCfg.getParameter<std::string>("optimisedPatsXmlFile");
  edm::LogImportant("PatternOptimizer") << " Writing optimized patterns to "<<fName << std::endl;
  XMLConfigWriter xmlWriter(omtfConfig);
  xmlWriter.writeGPs(goldenPatterns, fName);
}

const SimTrack* PatternOptimizer::findSimMuon(const edm::Event &event, const SimTrack * previous) {
  const SimTrack* result = 0;
  edm::Handle<edm::SimTrackContainer> simTks;
  event.getByLabel(edmCfg.getParameter<edm::InputTag>("g4SimTrackSrc"), simTks);

  for (std::vector<SimTrack>::const_iterator it=simTks->begin(); it< simTks->end(); it++) {
    const SimTrack& aTrack = *it;
    if ( !(aTrack.type() == 13 || aTrack.type() == -13) )
      continue;
    if(previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(), previous->momentum()) < 0.07)
      continue;
    if ( !result || aTrack.momentum().pt() > result->momentum().pt())
      result = &aTrack;
  }
  return result;
}

void PatternOptimizer::updateStat(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp, double delta, double norm) {
  //cout<<__FUNCTION__<<":"<<__LINE__<<" exptCandGp Pt "<<exptCandGp->key().thePt<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCandGp Pt "<<omtfCandGp->key().thePt<<" delta "<<delta<<" norm "<<norm<<std::endl;
  for(unsigned int iLayer = 0;  iLayer < omtfResult.getHitPdfBins().size(); iLayer++) {
    //updating statistic for the gp which should have fired
    int iBinExpt = exptResult.getHitPdfBins()[iLayer];
    if(iBinExpt != 0) {
      exptCandGp->updateStat(iLayer, exptResult.getRefLayer(), iBinExpt, whatExptVal, delta);
      exptCandGp->updateStat(iLayer, exptResult.getRefLayer(), iBinExpt, whatExptNorm, norm);
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
    }

    if(omtfResult.getHitPdfBins()[iLayer] != 0 ) {
      //updating statistic for the gp which found the candidate
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for omtf gp "<<omtfCandGp<<std::endl;
      omtfCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatOmtfVal, -delta);
      omtfCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatOmtfNorm, norm);
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for omtfCandGp: iLayer "<<iLayer<<" RefLayer "<<omtfCand.getRefLayer()<<" iBinOmtf "<<omtfResult.getHitPdfBins()[iLayer]<<" fired "<<omtfResult.isLayerFired(iLayer)<<std::endl;
    }
  }
}

void PatternOptimizer::updateStatVoter_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  double delta = 0;
  double norm = 1.;
  if(exptCandGp->key().thePt >= 51 && omtfCandGp->key().thePt < 51) {
    delta = 1.; //1./10.;
  }
  else if(exptCandGp->key().thePt <= 33 && omtfCandGp->key().thePt >= 51) {
    delta = 6./sqrt(exptCandGp->key().thePt);
    //norm = delta;
  }
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  /*if(exptCandGp->key().theCharge != omtfCandGp->key().theCharge) {
    delta = 2 * delta; //TODO what to do in this case????????
  }*/
  //delta = abs(delta);
  //delta *=delta;

  if(delta != 0) {
    updateStat(omtfCandGp, exptCandGp,  delta, norm);
  }
}

void PatternOptimizer::updateStatPtLogDiff_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  double delta = 0;
  if(exptCandGp->key().thePt > omtfCandGp->key().thePt) {
    delta = ( log(exptCandGp->key().thePt) - log(omtfCandGp->key().thePt) ) / log(omtfCandGp->key().thePt) ;
  }
  else {//exptCandGp->key().thePt <= omtfCandGp->key().thePt
    /*if(exptCandGp->key().thePt > 50)
      delta = 0;
    else*/
      delta = ( log(exptCandGp->key().thePt) - log(omtfCandGp->key().thePt) ) / log(exptCandGp->key().thePt) ;
  }
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  if(exptCandGp->key().theCharge != omtfCandGp->key().theCharge) {
    delta = 2 * delta; //TODO what to do in this case????????
  }
  //delta = abs(delta);
  delta *=delta;

  double norm = 1.;
  updateStat(omtfCandGp, exptCandGp,  delta, norm);
}

void PatternOptimizer::updateStatPtDiff2_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  double delta = 0;
  if(exptCandGp->key().thePt > omtfCandGp->key().thePt) {
    delta = ( exptCandGp->key().thePt - omtfCandGp->key().thePt);
    delta /= 10.;
  }
  else {//exptCandGp->key().thePt <= omtfCandGp->key().thePt
    delta = (omtfCandGp->key().thePt - exptCandGp->key().thePt);// / 2.; //watch out - the thePt is unsigned!!!!!!!
  }

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  if(omtfCandGp->key().thePt < 100 && exptCandGp->key().theCharge != omtfCandGp->key().theCharge) {
    delta = 2 * delta; //TODO what to do in this case????????
  }
  //double norm = 1./((double)exptCandGp->key().thePt * (double)exptCandGp->key().thePt);
  //double norm = exp(-1.0 * (double)exptCandGp->key().thePt);
  double norm = 1.;
  //delta = abs(delta);
  delta *=delta;
  //delta *= norm;

  if(omtfCandGp->key().thePt != exptCandGp->key().thePt)
    updateStat(omtfCandGp, exptCandGp,  delta, norm);
}


void PatternOptimizer::updateStatPtLogDiff_2(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  double delta = 0;
  if(exptCandGp->key().thePt > omtfCandGp->key().thePt) {
    delta = ( log(exptCandGp->key().thePt - omtfCandGp->key().thePt) ) / 1.5;
  }
  else {//exptCandGp->key().thePt <= omtfCandGp->key().thePt
    if(exptCandGp->key().thePt > 50)
      delta = 0;
    else
      delta = ( log(omtfCandGp->key().thePt) - log(exptCandGp->key().thePt) )  ;
  }
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  if(exptCandGp->key().theCharge != omtfCandGp->key().theCharge) {
    delta = 2 * delta; //TODO what to do in this case????????
  }
  //delta = abs(delta);
  //delta *=delta;
  double norm = 1./(double)exptCandGp->key().thePt;
  delta *= norm;
  updateStat(omtfCandGp, exptCandGp,  delta, norm);
}

void PatternOptimizer::updatePdfsMean_1(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) {
  for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
    double d = 0;
    if(gp->statisitics[iLayer][iRefLayer][iPdf][whatExptNorm] != 0) {
      d += gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptVal]/(double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm];
    }

    if(gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfNorm] != 0)
      d += gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfVal]/(double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm] ;


    d = d * learingRate;

    if(d > 0 && gp->key().thePt > 50 && gp->pdfAllRef[iLayer][iRefLayer][iPdf] == 0) {
      d = 0; //don't add additional non zero beans for the high pt
    }

    gp->pdfAllRef[iLayer][iRefLayer][iPdf] += d;

    if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] < 0)
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
    else if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] > omtfConfig->pdfMaxValue()) {
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = omtfConfig->pdfMaxValue();
    }


    if(d != 0) {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<" "<< gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" iBin "
          <<iPdf<<" pdfVal "<<gp->pdfAllRef[iLayer][iRefLayer][iPdf]
          <<" ExptVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][whatExptVal]
          <<" ExptNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][whatExptNorm]
          <<" OmtfVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfVal]
          <<" OmtfNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfNorm]
          <<" d "<<d<<std::endl;
    }
  }
}

void PatternOptimizer::updatePdfsMean_2(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) {
  for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
    double d = 0;
    double norm = (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm] + (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm];
    if(norm != 0)
      d += (gp->statisitics[iLayer][iRefLayer][iPdf][whatExptVal] + gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfVal])/norm;

    d = d * learingRate;

    if(d > 0 && gp->key().thePt > 50 && gp->pdfAllRef[iLayer][iRefLayer][iPdf] == 0) {
      d = 0; //don't add additional non zero beans for the high pt
    }

    gp->pdfAllRef[iLayer][iRefLayer][iPdf] += d;
    if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] < 0)
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
    else if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] > omtfConfig->pdfMaxValue()) {
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = omtfConfig->pdfMaxValue();
    }
    if(d != 0) {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<" "<< gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" iBin "
          <<iPdf<<" pdfVal "<<gp->pdfAllRef[iLayer][iRefLayer][iPdf]
          <<" ExptVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptVal]
          <<" ExptNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm]
          <<" OmtfVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfVal]
          <<" OmtfNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm]
          <<" d "<<d<<std::endl;
    }
  }
}

void PatternOptimizer::updatePdfsVoter_1(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) {
  for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
    double d = 0;
    double norm = (double)gp->statisitics[iLayer][iRefLayer][iPdf][whatExptNorm] + (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm];
    if(norm != 0)
      d += (gp->statisitics[iLayer][iRefLayer][iPdf][whatExptVal] + gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfVal]); // /norm;

    //d = d * learingRate;

    if(d > 0 && gp->key().thePt > 50 && gp->pdfAllRef[iLayer][iRefLayer][iPdf] == 0) {
      d = 0; //don't add additional non zero beans for the high pt
    }

    if(gp->statisitics[iLayer][iRefLayer][iPdf][whatExptNorm] > 5 || gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm] > 5) {
      if(d >= 2)
        d = 1.0;
      else if(d <= 2)
        d = -1.0;
    }
    else
      d = 0;

    gp->pdfAllRef[iLayer][iRefLayer][iPdf] += d;
    if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] < 0)
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
    else if(gp->pdfAllRef[iLayer][iRefLayer][iPdf] > omtfConfig->pdfMaxValue()) {
      gp->pdfAllRef[iLayer][iRefLayer][iPdf] = omtfConfig->pdfMaxValue();
    }
    if(d != 0) {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<" "<< gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" iBin "
          <<iPdf<<" pdfVal "<<gp->pdfAllRef[iLayer][iRefLayer][iPdf]
          <<" ExptVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptVal]
          <<" ExptNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm]
          <<" OmtfVal "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfVal]
          <<" OmtfNorm "<<gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm]
          <<" d "<<d<<std::endl;
    }
  }
}
