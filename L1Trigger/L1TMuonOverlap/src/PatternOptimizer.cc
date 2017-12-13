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

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"


PatternOptimizer::PatternOptimizer(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  edmCfg(edmCfg), omtfConfig(omtfConfig), goldenPatterns(gps), simMuon(0),
  //TODO set desire function here, see https://www.cprogramming.com/c++11/c++11-lambda-closures.html
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatCollectProb(omtfCandGp, exptCandGp); } ),
  updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatForAllGps(omtfCandGp, exptCandGp); } ),

  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatVoter_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtDiff_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtDiff2_1(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatPtLogDiff_2(omtfCandGp, exptCandGp); } ),
  //updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_1(gp, iLayer, iRefLayer, learingRate); })
  //updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsMean_2(gp, iLayer, iRefLayer, learingRate); })
  updatePdfsFunc([this] (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) { updatePdfsVoter_1(gp, iLayer, iRefLayer, learingRate); })
{
  //modifyPatterns(); //TODO remove if not needed!!!!!!!!!!!!!!!!
  simMuPt =  new TH1I("simMuPt", "simMuPt", goldenPatterns.size(), -0.5, goldenPatterns.size()-0.5);
  simMuFoundByOmtfPt =  new TH1I("simMuFoundByOmtfPt", "simMuFoundByOmtfPt", goldenPatterns.size(), -0.5, goldenPatterns.size()-0.5);
}

void PatternOptimizer::modifyPatterns() {
  cout<<__FUNCTION__<<": "<<__LINE__<<" called!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          if(gp->key().thePt > 30 && gp->pdfAllRef[iLayer][iRefLayer][iPdf] < 0.01) //suppressing tails of the distributions
            gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
        }
      }
    }
  }
}

PatternOptimizer::~PatternOptimizer() {
}

void PatternOptimizer::observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
    const std::vector<AlgoMuon>& algoCandidates,
    std::vector<AlgoMuon>& gbCandidates,
    const std::vector<l1t::RegionalMuonCand> & candMuons) {

  unsigned int procIndx = omtfConfig->getProcIndx(iProcessor, mtfType);

/*
  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;
  int patNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(patNum).get(); // expected pattern
*/

  //bool found = false;
  for(auto& gbCandidate : gbCandidates) {
    int iRefHit = gbCandidate.getRefHitNumber();
    if(gbCandidate.getGoldenPatern() != 0 && gbCandidate.getGoldenPatern()->getResults()[procIndx][iRefHit].getFiredLayerCnt() > omtfResult.getFiredLayerCnt() ) {
      //cout<<__FUNCTION__<<":"<<__LINE__<<" gbCandidate "<<gbCandidate<<" "<<std::endl;
      omtfCand = gbCandidate;
      omtfResult = gbCandidate.getGoldenPatern()->getResults()[procIndx][iRefHit];
      //exptResult = exptCandGp->getResults()[procIndx][iRefHit];
      candProcIndx = procIndx;
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
        if(omtfResult.getFiredLayerCnt() == gp->getResults()[procIndx][iRefHit].getFiredLayerCnt() )
        {
          cout<<gp->key()<<std::endl<<gp->getResults()[procIndx][iRefHit]<<std::endl;
        }
      }
      std::cout<<std::endl;
    }
  }*/
}

void PatternOptimizer::observeEventBegin(const edm::Event& iEvent) {
  omtfCand = AlgoMuon();
  candProcIndx = 0xffff;
  omtfResult =  GoldenPatternResult();
  //exptResult =  GoldenPatternResult();

  simMuon = findSimMuon(iEvent);
  //cout<<__FUNCTION__<<":"<<__LINE__<<" simMuon "<<simMuon->momentum().pt()<<" "<<std::endl;
}

void PatternOptimizer::observeEventEnd(const edm::Event& iEvent) {
  if(simMuon == 0 || omtfCand.getGoldenPatern() == 0)//no sim muon or empty candidate
    return;

  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCand "<<omtfCand<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfResult "<<std::endl<<omtfResult<<std::endl;

  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;

  GoldenPatternWithStat* omtfCandGp = static_cast<GoldenPatternWithStat*>(omtfCand.getGoldenPatern());

  exptPatNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(exptPatNum).get(); // expected pattern

  int iRefHit = omtfCand.getRefHitNumber();
  //GoldenPatternResult omtfResult = omtfCandGp->getResults()[candProcIndx][iRefHit];
  exptResult = exptCandGp->getResults()[candProcIndx][iRefHit];

  simMuFoundByOmtfPt->Fill(exptCandGp->key().theNumber); //TODO add weight of the muons pt spectrum

  updateStatFunc(omtfCandGp, exptCandGp);

  ///debug printout
  if( omtfCandGp->key().thePt > 40 && exptCandGp->key().thePt <= 15 )
  {
    std::cout<<iEvent.id()<<std::endl;
    std::cout<<" ptSim "<<ptSim<<" chargeSim "<<chargeSim<<" patNum "<<exptPatNum<<std::endl;
    std::cout<<"exptCandGp "<<exptCandGp->key()<<std::endl;
    std::cout<<"exptResult "<<std::endl<<exptResult<<std::endl;

    std::cout<<"omtfCandGp "<<omtfCandGp->key()<<std::endl;
    std::cout<<"omtfResult "<<std::endl<<omtfResult<<std::endl;
    /*std::cout<<"other gps results"<<endl;
    for(auto& gp : goldenPatterns) {
      if(omtfResult.getFiredLayerCnt() == gp->getResults()[candProcIndx][iRefHit].getFiredLayerCnt() )
      {
        cout<<gp->key()<<std::endl<<gp->getResults()[candProcIndx][iRefHit]<<std::endl;
      }
    }*/
    std::cout<<std::endl;
  }

}

void PatternOptimizer::endJob() {
  //savePatternsInRoot("orginalPatterns.root");

  //calulateProb();

/*  double learingRate = 0.1;
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        updatePdfsFunc(gp.get(), iLayer, iRefLayer, learingRate);
      }
    }
  }*/

  //"optimisedPats.xml";
  std::string fName = edmCfg.getParameter<std::string>("optimisedPatsXmlFile");
  edm::LogImportant("PatternOptimizer") << " Writing optimized patterns to "<<fName << std::endl;
  XMLConfigWriter xmlWriter(omtfConfig, true, false);
  xmlWriter.writeGPs(goldenPatterns, fName);

  fName.replace(fName.find('.'), fName.length(), ".root");
  savePatternsInRoot(fName);
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


void PatternOptimizer::updateStatCollectProb(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
/*  for(unsigned int iLayer = 0;  iLayer < exptResult.getHitPdfBins().size(); iLayer++) {
    //updating statistic for the gp which should have fired
    int iBinExpt = exptResult.getHitPdfBins()[iLayer];
    if(iBinExpt != 0) {
      exptCandGp->updateStat(iLayer, exptResult.getRefLayer(), iBinExpt, whatExptVal, 1); //TODO in principle the events shoud be weighted by the muon pt probability
      exptCandGp->updateStat(iLayer, exptResult.getRefLayer(), iBinExpt, whatExptNorm, 1);
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
    }
  }*/

  for(unsigned int iRefHit = 0; iRefHit < exptCandGp->getResults()[candProcIndx].size(); iRefHit++) {
    GoldenPatternResult result = exptCandGp->getResults()[candProcIndx][iRefHit];

    for(unsigned int iLayer = 0;  iLayer < result.getHitPdfBins().size(); iLayer++) {
      //updating statistic for the gp which should have fired
      int iBinExpt = result.getHitPdfBins()[iLayer];
      if(iBinExpt != 0) {
        exptCandGp->updateStat(iLayer, result.getRefLayer(), iBinExpt, whatExptVal, 1); //TODO in principle the events should be weighted by the muon pt probability
        exptCandGp->updateStat(iLayer, result.getRefLayer(), iBinExpt, whatExptNorm, 1);
        //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
      }
    }
  }
}

void PatternOptimizer::calulateProb() {
  for(auto& gp : goldenPatterns) {
    cout<<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<" Calculating P(x | C_k) "<<std::endl;
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
        if(iLayer == refLayerLogicNumber) //skip this as here we keep the P(C_k),
          continue;

        double pdfNorm = 0;
        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          pdfNorm += (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm];
        }

        cout<<__FUNCTION__<<":"<<__LINE__<<" : iLayer "<<iLayer<<" RefLayer "<<iRefLayer<<" pdfNorm "<<pdfNorm<<std::endl;

        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          if(pdfNorm > 50) {//50 is to reject the one with very small statistics
            double prob = (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm] / pdfNorm;
            gp->pdfAllRef[iLayer][iRefLayer][iPdf] = prob;
          }
          else
            gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
        }
      }
    }
  }

  //in pdf value for the ref layer we keep the "class probability" P(C_k), where class is the pt-sign (i.e. golde pattern)
  cout<<__FUNCTION__<<": "<<__LINE__<<" Calculating P(C_k) "<<std::endl;
  unsigned int iPdf = omtfConfig->nPdfBins()/2;// <<(omtfConfig->nPdfAddrBits()-1);
  for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns[0]->getPdf()[0].size(); ++iRefLayer) {
    double norm = 0;
    unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
    for(auto& gp : goldenPatterns) {
      norm += (double)gp->statisitics[refLayerLogicNumber][iRefLayer][iPdf][ whatExptNorm];
    }

    cout<<__FUNCTION__<<":"<<__LINE__<<" RefLayer "<<iRefLayer<<" norm "<<norm<<std::endl;
    for(auto& gp : goldenPatterns) {
      if(norm > 50)
        gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf] = (double)gp->statisitics[refLayerLogicNumber][iRefLayer][iPdf][ whatExptNorm] / norm;
      else
        gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf] = 0;
    }
  }
}


void PatternOptimizer::updateStatForAllGps(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  unsigned int iRefHit = omtfCand.getRefHitNumber();
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;
  for(auto& itGP: goldenPatterns) {
    if(itGP->key().thePt == 0 )
      continue;
    auto& result = itGP->getResults()[candProcIndx][iRefHit];
    if(result.getFiredLayerCnt() < omtfResult.getFiredLayerCnt())
      continue;
    /*//this loop will be needed for training
     * for(unsigned int iLayer = 0;  iLayer < result.getHitPdfBins().size(); iLayer++) {
      int iBin = result.getHitPdfBins()[iLayer];
      if(iBin != 0) {
        itGP->updateStat(iLayer, result.getRefLayer(), iBin, whatExptVal, 1);
        itGP->updateStat(iLayer, result.getRefLayer(), iBin, whatExptNorm, 1);
        //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
      }
    }*/
    double gpProbability = result.getGpProbability1();
    //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<itGP->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<" FiredLayerCnt "<<result.getFiredLayerCnt()<<std::endl;
    itGP->gpProbabilityStat[exptPatNum].Fill(gpProbability);
  }
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



void PatternOptimizer::savePatternsInRoot(std::string rootFileName) {
  gStyle->SetOptStat(111111);
  TFile outfile(rootFileName.c_str(), "RECREATE");
  cout<<__FUNCTION__<<": "<<__LINE__<<" out fileName "<<rootFileName<<" outfile->GetName() "<<outfile.GetName()<<endl;

  outfile.cd();
  simMuFoundByOmtfPt->Write();

  outfile.mkdir("patternsPdfs")->cd();
  ostringstream ostrName;
  ostringstream ostrTtle;
  for(auto& gp : goldenPatterns) {
    OMTFConfiguration::PatternPt patternPt = omtfConfig->getPatternPtRange(gp->key().theNumber);
    if(gp->key().thePt == 0)
      continue;
    //cout<<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<std::endl;
    ostrName.str("");
    ostrName<<"PatNum_"<<gp->key().theNumber;
    ostrTtle.str("");
    ostrTtle<<"PatNum_"<<gp->key().theNumber<<"_ptCode_"<<gp->key().thePt<<"_Pt_"<<patternPt.ptFrom<<"_"<<patternPt.ptTo<<"_GeV";
    TCanvas* canvas = new TCanvas(ostrName.str().c_str(), ostrTtle.str().c_str(), 1200, 1000);
    canvas->Divide(gp->getPdf().size(), gp->getPdf()[0].size(), 0, 0);
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        canvas->cd(1 + iLayer + iRefLayer * gp->getPdf().size());
        //unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
        ostrName.str("");
        ostrName<<"PatNum_"<<gp->key().theNumber<<"_refLayer_"<<iRefLayer<<"_Layer_"<<iLayer;
        ostrTtle.str("");
        ostrTtle<<"PatNum "<<gp->key().theNumber<<" ptCode "<<gp->key().thePt<<" refLayer "<<iRefLayer<<" Layer "<<iLayer<<" meanDistPhi "<<gp->meanDistPhi[iLayer][iRefLayer][0]; //"_Pt_"<<patternPt.ptFrom<<"_"<<patternPt.ptTo<<"_GeV
        //cout<<__FUNCTION__<<": "<<__LINE__<<" creating hist "<<ostrTtle.str()<<std::endl;
        TH1F* hist = new TH1F(ostrName.str().c_str(), ostrTtle.str().c_str(), omtfConfig->nPdfBins(), -0.5, omtfConfig->nPdfBins()-0.5);
        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          hist->Fill(iPdf, gp->pdfAllRef[iLayer][iRefLayer][iPdf]);
        }
        hist->Write();
        hist->Draw("hist");
      }
    }
    canvas->Write();
    delete canvas;
  }

  outfile.mkdir("patternsPdfSumStat")->cd();
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0 )
      continue;
    for(auto& ptBin : gp->gpProbabilityStat) {
      if(ptBin.GetNbinsX() > 5)
        ptBin.Write();
    }
  }
  outfile.Close();
}
