/*
 * PatternOptimizer.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PatternOptimizer.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"

#include "Math/VectorUtil.h"

PatternOptimizer::PatternOptimizer(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  edmCfg(edmCfg), omtfConfig(omtfConfig), goldenPatterns(gps), simMuon(0),
  //TODPO set desire function here, see https://www.cprogramming.com/c++11/c++11-lambda-closures.html
  updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtLogDiff_1(omtfCandGp, exptCandGp); } ),
  //updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_1(gp, iLayer, iRefLayer, learingRate); })
  updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_2(gp, iLayer, iRefLayer, learingRate); })
{

}

PatternOptimizer::~PatternOptimizer() {
}

void PatternOptimizer::observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
    const std::vector<AlgoMuon>& algoCandidates,
    std::vector<AlgoMuon>& gbCandidates,
    const std::vector<l1t::RegionalMuonCand> & candMuons) {

  for(auto& gbCandidate : gbCandidates) {
    int refHitNum = gbCandidate.getRefHitNumber();
    if(gbCandidate.getGoldenPatern() != 0 && gbCandidate.getGoldenPatern()->getResults().at(refHitNum).getFiredLayerCnt() > omtfResult.getFiredLayerCnt() ) {
      //cout<<__FUNCTION__<<":"<<__LINE__<<" gbCandidate "<<gbCandidate<<" "<<std::endl;
      omtfCand = gbCandidate;
      omtfResult = gbCandidate.getGoldenPatern()->getResults().at(refHitNum);
    }
  }
}

void PatternOptimizer::observeEventBegin(const edm::Event& iEvent) {
  omtfCand = AlgoMuon();
  omtfResult =  GoldenPatternResult();

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
  //cout<<__FUNCTION__<<":"<<__LINE__<<" ptSim "<<ptSim<<" patNum "<<patNum<<" exptCandGp "<<exptCandGp<<std::endl;

  updateStatFunc(omtfCandGp, exptCandGp);

}

void PatternOptimizer::endJob() {
  double learingRate = 1.;
  for(auto& gp : goldenPatterns) {
    //gp->updatePdfs(learingRate);

    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        updatePdfsFunc(gp.get(), iLayer, iRefLayer, learingRate);
      }
    }
  }

  std::string fName = "optimisedPats.xml";
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
  delta = abs(delta);
  //delta *=delta;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" exptCandGp Pt "<<exptCandGp->key().thePt<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCandGp Pt "<<omtfCandGp->key().thePt<<" delta "<<delta<<std::endl;

  for(unsigned int iLayer = 0;  iLayer < omtfResult.getHitPdfBins().size(); iLayer++) {
    if(omtfResult.isLayerFired(iLayer) ) {
      //updating statistic for the gp which should have fired
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for sim gp "<<std::endl;
      exptCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatExptVal, delta);
      exptCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatExptNorm, 1);

      //updating statistic for the gp which found the candidate
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for omtf gp "<<omtfCandGp<<std::endl;
      omtfCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatOmtfVal, -delta);
      omtfCandGp->updateStat(iLayer, omtfCand.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], whatOmtfNorm, 1);
    }
  }
}

void PatternOptimizer::updatePdfsMean_1(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) {
  for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
    double d = 0;
    if(gp->statisitics[iLayer][iRefLayer][iPdf][whatExptNorm] != 0)
      d += gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptVal]/(double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm];

    if(gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfNorm] != 0)
      d += gp->statisitics[iLayer][iRefLayer][iPdf][whatOmtfVal]/(double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatOmtfNorm] ;

    d = d * learingRate;
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
