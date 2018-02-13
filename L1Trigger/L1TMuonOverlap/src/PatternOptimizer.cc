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

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include "Math/VectorUtil.h"

#include <boost/range/adaptor/reversed.hpp>

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"


PatternOptimizer::PatternOptimizer(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  edmCfg(edmCfg), omtfConfig(omtfConfig), goldenPatterns(gps), simMuon(0),
  //TODO set desire function here, see https://www.cprogramming.com/c++11/c++11-lambda-closures.html
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatCollectProb(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatForAllGps(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { calculateThresholds(omtfCandGp, exptCandGp); } ),
  //updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { tuneClassProb(omtfCandGp, exptCandGp); } ),
  updateStatFunc([this] (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) { updateStatCloseResults(omtfCandGp, exptCandGp); } ),

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

  //ptRangeFrom = edmCfg.getParameter<double>("ptRangeFrom");
  //ptRangeFrom = edmCfg.getParameter<double>("ptRangeTo");

  selectedPatNum = edmCfg.getParameter<unsigned int>("selectedPatNum"); //TODO

  ptCut = goldenPatterns[selectedPatNum]->key().thePt;

  currnetPtBatchPatNum = selectedPatNum;

  //printPatterns();

  initRateWeights();
}

void PatternOptimizer::modifyPatterns() {
  cout<<__FUNCTION__<<": "<<__LINE__<<" called!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;
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

void PatternOptimizer::printPatterns() {
  cout<<__FUNCTION__<<": "<<__LINE__<<" called!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  for(int patNum = goldenPatterns.size() -1; patNum >= 0; patNum--  ) {
    double pt = omtfConfig->getPatternPtRange(patNum).ptFrom;
    if(pt > 0) {
      cout<<"cmsRun runThresholdCalc.py "<<patNum<<" "<<(patNum+1)<<" _"<<RPCConst::iptFromPt(pt)<<"_";
      if(goldenPatterns[patNum]->key().theCharge == -1)
        cout<<"m_";
      else
        cout<<"p_";

      cout<<" > out"<<patNum<<".txt"<<std::endl;
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
  //cout<<__FUNCTION__<<":"<<__LINE__<<" evevt "<<iEvent.id().event()<<" simMuon pt "<<simMuon->momentum().pt()<<" GeV "<<std::endl;
}

void PatternOptimizer::observeEventEnd(const edm::Event& iEvent) {
  if(simMuon == 0 || omtfCand.getGoldenPatern() == 0)//no sim muon or empty candidate
    return;

  //cout<<__FUNCTION__<<":"<<__LINE__<<" event "<<iEvent.id().event()<<endl;
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

  //TODO select the function to be executed
  updatePdfCloseResults();

  //calulateProb();

/*  double learingRate = 0.1;
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        updatePdfsFunc(gp.get(), iLayer, iRefLayer, learingRate);
      }
    }
  }*/

/*
  double step = edmCfg.getParameter<double>("step");
  modifyPatterns1(step); //TODO remove!!!!!!!!!!!!!!!!!!!!!!!!
*/

  std::string fName = edmCfg.getParameter<std::string>("optimisedPatsXmlFile");
  edm::LogImportant("PatternOptimizer") << " Writing optimized patterns to "<<fName << std::endl;
  XMLConfigWriter xmlWriter(omtfConfig, true, false);
  xmlWriter.writeGPs(goldenPatterns, fName);

  //calculateThresholds
/*  for(int i = 1; i <= 1; i++) {
    //TODO select what to do
    double targetEff = 0.795 + 0.0005 * i;
    calculateThresholds(targetEff);

    double targetEff = edmCfg.getParameter<double>("targetEff"); //0.01 * i; //targetEffLost
    tuneClassProb(targetEff); //connot be done in the loop, because the pdf's are updated each time, to make it good the original gp should be copied

    std::ostringstream fName;// = edmCfg.getParameter<std::string>("optimisedPatsXmlFile");

    //fName<<"optimisedPats_thr_0"<<setw(4)<<setfill('0')<<(targetEff * 10000)<<".xml";
    fName<<"optimisedPats_"<<selectedPatNum<<".xml";

    edm::LogImportant("PatternOptimizer") << " Writing  patterns with thresholds to "<<fName.str() << std::endl;
    //XMLConfigWriter xmlWriter(omtfConfig, true, false);
    xmlWriter.writeGPs(goldenPatterns, fName.str());
  }*/

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
        exptCandGp->updateStat(iLayer, result.getRefLayer(), iBinExpt, whatExptNorm, rateWeights[exptPatNum]);
        //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
      }
    }
  }
}

void PatternOptimizer::calulateProb() {
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;
    cout<<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<" Calculating P(x | C_k) "<<std::endl;
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
        if(iLayer == refLayerLogicNumber) //skip this as here we keep the P(C_k),
          continue;

        double pdfCnt = 0;
        double pdfNorm = 0;
        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          pdfCnt  += (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptVal];
          pdfNorm += (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm];
        }

        cout<<__FUNCTION__<<":"<<__LINE__<<" : iLayer "<<iLayer<<" RefLayer "<<iRefLayer<<" pdfNorm "<<pdfNorm<<std::endl;

        for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          if(pdfCnt > 50) {//50 is to reject the one with very small statistics
            double prob = (double)gp->statisitics[iLayer][iRefLayer][iPdf][ whatExptNorm] / pdfNorm;
            gp->pdfAllRef[iLayer][iRefLayer][iPdf] = prob;
          }
          else
            gp->pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
        }
      }
    }
  }

  //in pdf[32] value for the ref layer we keep the "class probability" P(C_k), where class is the pt-sign (i.e. golde pattern)
  cout<<__FUNCTION__<<": "<<__LINE__<<" Calculating P(C_k) "<<std::endl;
  unsigned int iPdf = omtfConfig->nPdfBins()/2;// <<(omtfConfig->nPdfAddrBits()-1);
  for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns[0]->getPdf()[0].size(); ++iRefLayer) {
    double cnt = 0;
    double norm = 0;
    unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
    for(auto& gp : goldenPatterns) {
      if(gp->key().thePt == 0)
        continue;
      cnt  += (double)gp->statisitics[refLayerLogicNumber][iRefLayer][iPdf][ whatExptVal];
      norm += (double)gp->statisitics[refLayerLogicNumber][iRefLayer][iPdf][ whatExptNorm];
    }

    cout<<__FUNCTION__<<":"<<__LINE__<<" RefLayer "<<iRefLayer<<" norm "<<norm<<std::endl;
    for(auto& gp : goldenPatterns) {
      if(gp->key().thePt == 0)
        continue;
      if(cnt > 50)
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
      itGP->gpProbabilityStat[result.getRefLayer()]->Fill(0);
    else {
      double gpProbability = result.getGpProbability1();
      //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<itGP->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<" FiredLayerCnt "<<result.getFiredLayerCnt()<<std::endl;
      itGP->gpProbabilityStat[result.getRefLayer()]->Fill(gpProbability);

      /*//this loop will be needed for training
       * for(unsigned int iLayer = 0;  iLayer < result.getHitPdfBins().size(); iLayer++) {
      int iBin = result.getHitPdfBins()[iLayer];
      if(iBin != 0) {
        itGP->updateStat(iLayer, result.getRefLayer(), iBin, whatExptVal, 1);
        itGP->updateStat(iLayer, result.getRefLayer(), iBin, whatExptNorm, 1);
        //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for exptCandGp: iLayer "<<iLayer<<" RefLayer "<<exptResult.getRefLayer()<<" iBinExpt "<<iBinExpt<<" fired "<<exptResult.isLayerFired(iLayer)<<std::endl;
      }
    }*/
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////

void PatternOptimizer::updateStatCloseResults(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  unsigned int iRefHit = omtfCand.getRefHitNumber();
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;

  if(omtfCandGp->key().theCharge != exptCandGp->key().theCharge) {
    /*cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCandGp->key().theCharge != exptCandGp->key().theCharge"<<endl;
    cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
    cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;
    cout<<__FUNCTION__<<":"<<__LINE__<<" --------------------------------"<<std::endl;*/
    return;
  }

  int omtfCandPtCode = patternPtCodes[omtfCandGp->key().theNumber];
  int exptCandPtCode = patternPtCodes[exptCandGp->key().theNumber];

  GoldenPatternWithStat* secondBestGp = 0;
  for(auto& itGP: goldenPatterns) {
    if(itGP->key().thePt == 0 || omtfCandGp->key().theCharge != itGP->key().theCharge)
      continue;
    auto& result = itGP->getResults()[candProcIndx][iRefHit];
    if(result.getFiredLayerCnt() < omtfResult.getFiredLayerCnt()) {
      //itGP->gpProbabilityStat.Fill(0);
    }
    else {
      //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<itGP->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<" FiredLayerCnt "<<result.getFiredLayerCnt()<<std::endl;
      //itGP->gpProbabilityStat.Fill(result.getPdfSum());

      //finding second highest results, for which PdfSum differs from the omtfCandGp by 1 ,
      if(itGP.get() != omtfCandGp) {
        double delta = 1;
        if(itGP->key().theNumber > omtfCandGp->key().theNumber)
          delta = 0;

        int thisGpPtCode = patternPtCodes[itGP->key().theNumber];
        if(omtfResult.getPdfSum() - result.getPdfSum() == delta) {
          if(secondBestGp != 0) {
            if(abs(patternPtCodes[secondBestGp->key().theNumber] - exptCandPtCode) > abs(thisGpPtCode - exptCandPtCode) ) {
              /*cout<<__FUNCTION__<<":"<<__LINE__<<" changing secondBestGp from "
                  <<secondBestGp->key()<<"\n"<<secondBestGp->getResults()[candProcIndx][iRefHit]<<"\nto "
                  <<itGP->key()<<result<<"\n omtfResult "
                  <<omtfCandGp->key()<<"\n"<<omtfResult<<"\n"
                  <<"exptCandGp"<<exptCandGp->key()<<endl;*/
              secondBestGp = itGP.get();
              // selecting secondBestGp that is closer to the exptCandGp
            }
          }
          else {
            secondBestGp = itGP.get();
          }

          //updating stat for any gp if the omtfCand is closer to the exptCand then this gp
          if( abs(omtfCandPtCode - exptCandPtCode) < abs(thisGpPtCode - exptCandPtCode) ) {
            for(unsigned int iLayer = 0;  iLayer < omtfResult.getHitPdfBins().size(); iLayer++) {
              itGP->updateStat(iLayer, result.getRefLayer(), result.getHitPdfBins()[iLayer], goodSmaller, 1);
            }
          }
        }
      }
    }
  }

  if(secondBestGp) {
    auto& secondBestResult = secondBestGp->getResults()[candProcIndx][iRefHit];
    int secondBestGpPtCode = patternPtCodes[secondBestGp->key().theNumber];

    if( abs(omtfCandPtCode - exptCandPtCode) < abs(secondBestGpPtCode - exptCandPtCode) ) {
      for(unsigned int iLayer = 0;  iLayer < omtfResult.getHitPdfBins().size(); iLayer++) {
        omtfCandGp->updateStat(iLayer, omtfResult.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], goodBigger, 1);
        //secondBestGp->updateStat(iLayer, secondBestResult.getRefLayer(), secondBestResult.getHitPdfBins()[iLayer], goodSmaller, 1);
      }
    }
    else if( (omtfCandPtCode - exptCandPtCode) > 2 && abs(omtfCandPtCode - exptCandPtCode) > abs(secondBestGpPtCode - exptCandPtCode) ) {
      for(unsigned int iLayer = 0;  iLayer < omtfResult.getHitPdfBins().size(); iLayer++) {
        //cout<<__FUNCTION__<<":"<<__LINE__<<" badBigger "<<iLayer<<" "<<omtfResult.getRefLayer()<<" "<<omtfResult.getHitPdfBins()[iLayer]<<endl;
        //cout<<__FUNCTION__<<":"<<__LINE__<<" badSmaller "<<iLayer<<" "<<secondBestResult.getRefLayer()<<" "<<secondBestResult.getHitPdfBins()[iLayer]<<endl;
        omtfCandGp->updateStat(iLayer, omtfResult.getRefLayer(), omtfResult.getHitPdfBins()[iLayer], badBigger, 1);
        secondBestGp->updateStat(iLayer, secondBestResult.getRefLayer(), secondBestResult.getHitPdfBins()[iLayer], badSmaller, 1);
      }
      /*cout<<__FUNCTION__<<":"<<__LINE__<<" omtf cand id bad. secondBestGp "
          <<secondBestGp->key()<<"\n"<<secondBestGp->getResults()[candProcIndx][iRefHit]
          <<"\n omtfResult "<<omtfCandGp->key()<<"\n"<<omtfResult<<"\n"
          <<"exptCandGp"<<exptCandGp->key()<<endl<<"---------------------------------------------"<<endl;*/
    }
    else {//equall
      //todo
    }

  }
}


void PatternOptimizer::updatePdfCloseResults() {
  cout<<__FUNCTION__<<":"<<__LINE__<<endl;
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;

    for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      unsigned int selectedBiggerLayer = 0;
      unsigned int selectedBiggerBin = 0;
      int maxDiffBigger = 0;

      unsigned int selectedSmallerLayer = 0;
      unsigned int selectedSmallerBin = 0;
      int maxDiffSmaller = 0;

      unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
       /* if(refLayerLogicNumber == iLayer)
          continue; //not taking into account the class probability*/

        for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
          int diffBigger = gp->statisitics[iLayer][iRefLayer][iPdf][badBigger] - gp->statisitics[iLayer][iRefLayer][iPdf][goodBigger];
          if(diffBigger > maxDiffBigger) {
            selectedBiggerLayer = iLayer;
            selectedBiggerBin = iPdf;
            maxDiffBigger = diffBigger;
          }

          int diffSmaller = gp->statisitics[iLayer][iRefLayer][iPdf][badSmaller] - gp->statisitics[iLayer][iRefLayer][iPdf][goodSmaller];
          if(diffSmaller > maxDiffSmaller) {
            selectedSmallerLayer = iLayer;
            selectedSmallerBin = iPdf;
            maxDiffSmaller = diffSmaller;
          }
        }
      }

      if(selectedBiggerBin != 0 && maxDiffBigger > 5) {
        gp->pdfAllRef[selectedBiggerLayer][iRefLayer][selectedBiggerBin]--;
        cout<<__FUNCTION__<<":"<<__LINE__<<" decreasing the pdf "<<gp->key()<<" iRefLayer "<<iRefLayer<<" layer "<<selectedBiggerLayer<<" bin "<<selectedBiggerBin
            <<" badBigger "<<gp->statisitics[selectedBiggerLayer][iRefLayer][selectedBiggerBin][badBigger]<<" goodBigger "<<gp->statisitics[selectedBiggerLayer][iRefLayer][selectedBiggerBin][goodBigger]
            <<" maxDiffBigger "<<maxDiffBigger<<endl;
      }

      if(selectedSmallerBin != 0 && maxDiffSmaller > 5) {
        gp->pdfAllRef[selectedSmallerLayer][iRefLayer][selectedSmallerBin]++;
        cout<<__FUNCTION__<<":"<<__LINE__<<" increasing the pdf "<<gp->key()<<" iRefLayer "<<iRefLayer<<" layer "<<selectedSmallerLayer<<" bin "<<selectedSmallerBin
            <<" badSmaller "<<gp->statisitics[selectedSmallerLayer][iRefLayer][selectedSmallerBin][badSmaller]<<" goodSmaller "<<gp->statisitics[selectedSmallerLayer][iRefLayer][selectedSmallerBin][goodSmaller]
            <<" maxDiffSmaller "<<maxDiffSmaller<<endl;
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////////////


void PatternOptimizer::calculateThresholds(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  unsigned int iRefHit = omtfCand.getRefHitNumber();
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;

  double gpProbability = exptResult.getGpProbability2(); //TODO chose  getGpProbability1 or getGpProbability2
  exptCandGp->gpProbabilityStat[exptResult.getRefLayer()]->Fill(gpProbability);
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<std::endl;
  return;

  ///////////////other version
  ///////////////////////////
  ///////////////////////////
  if(selectedPatNum != exptPatNum) {
    return;
  }
  if(omtfCandGp->key().thePt > exptCandGp->key().thePt) {
    exptCandGp->gpProbabilityStat[exptResult.getRefLayer()]->Fill(1.1); //filling overflow bin in this case
  }
  else if(omtfCandGp == exptCandGp) {// for the exptCandGp the threshold should be 0 in this iteration, so the ghosBuster should choose it, if the gp with the higher pt wsa not chose
    //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<std::endl;
    exptCandGp->gpProbabilityStat[exptResult.getRefLayer()]->Fill(gpProbability);

    //debug
    /*if( gpProbability <= 0.002 ) {
      std::cout<<"omtfResult "<<omtfCandGp->key()<<std::endl;
      std::cout<<std::endl<<omtfResult<<std::endl;
      //int refHitNum = omtfCand.getRefHitNumber();
      std::cout<<"other gps results"<<endl;
      for(auto& gp : goldenPatterns) {
        if(omtfResult.getFiredLayerCnt() == gp->getResults()[candProcIndx][iRefHit].getFiredLayerCnt() )
        {
          cout<<gp->key()<<std::endl<<gp->getResults()[candProcIndx][iRefHit]<<std::endl;
        }
      }
      std::cout<<std::endl;
    }*/

  }
  else {//it can happened if the exptCandGp has lower FiredLayerCnt then the omtfCandGp
    exptCandGp->gpProbabilityStat[exptResult.getRefLayer()]->Fill(-0.1); //filling underflow bin in this case

/*    auto& result = exptCandGp->getResults()[candProcIndx][iRefHit];
    double gpProbability = exptCandGp->getResults()[candProcIndx][iRefHit].getGpProbability1();
    cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" iRefHit "<<iRefHit<<" gpProbability "<<gpProbability<<" FiredLayerCnt "<<result.getFiredLayerCnt()<<std::endl;*/
  }
}

/*
void PatternOptimizer::calculateThresholds(double targetEff) {
  cout<<__FUNCTION__<<":"<<__LINE__<<" targetEff "<<targetEff<<std::endl;
/*  TFile outfile("optimisedPats_2.root", "READ"); //FIXME the file name

  ostringstream ostrName;
  ostringstream ostrTitle;
  for(unsigned int iPat = 0; iPat < this->myOmtfConfig->nGoldenPatterns(); iPat++) {
    ostrName.str("");
    ostrTitle.str("");
    OMTFConfiguration::PatternPt patternPt = this->myOmtfConfig->getPatternPtRange(iPat);
    ostrName<<"gpProbabilityStat_GP_"<<key().theNumber<<"_ptBinNum_"<<iPat;//<<"_muPtFrom_"<<patternPt.ptFrom<<"_GeV";
    ostrTitle<<"gpProbabilityStat_GP_"<<key().theNumber<<"_ptBinNum_"<<iPat<<"_muPtFrom_"<<patternPt.ptFrom<<"_GeV";
    if(patternPt.ptFrom > 0)
      gpProbabilityStat.emplace_back(TH1I(ostrName.str().c_str(), ostrTitle.str().c_str(), 100, 0., 1.)); //TODO find proper range
    else
      gpProbabilityStat.emplace_back(TH1I(ostrName.str().c_str(), ostrTitle.str().c_str(), 1, 0., 1.)); //to save some memory, for empty patterns just "empty" hits
  }
  *

  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;

  int takenMuCnt = 0;
  for(int iBin = goldenPatterns[currnetPtBatchPatNum]->gpProbabilityStat.GetNbinsX() + 1; iBin >= 0; iBin--) {//we include overflow bin, because it contains the muons found by the gps with higher pt
    takenMuCnt += goldenPatterns[currnetPtBatchPatNum]->gpProbabilityStat.GetBinContent(iBin);
    if(takenMuCnt >= targetEff * (double)(simMuFoundByOmtfPt->GetBinContent(currnetPtBatchPatNum + 1)) ) {
      float threshold = goldenPatterns[currnetPtBatchPatNum]->gpProbabilityStat.GetXaxis()->GetBinLowEdge(iBin);
      goldenPatterns[currnetPtBatchPatNum]->setThreshold(0, threshold);
      cout<<__FUNCTION__<<":"<<__LINE__<<" currnetPtBatchPatNum "<<currnetPtBatchPatNum<<" takenMuCnt "<<takenMuCnt
          <<" simMuFoundByOmtfPt "<<simMuFoundByOmtfPt->GetBinContent(currnetPtBatchPatNum + 1)
          <<" "<<goldenPatterns[currnetPtBatchPatNum]->gpProbabilityStat.Integral(0, 100000)
          <<" found threshold "<<threshold<<std::endl;
      break;
    }
  }

}
*/

void PatternOptimizer::calculateThresholds(double targetEff) {
  cout<<__FUNCTION__<<":"<<__LINE__<<" targetEff "<<targetEff<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;

  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;
    for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      int takenMuCnt = 0;
      int allMuCnt = gp->gpProbabilityStat[iRefLayer]->Integral(0, 100000);

      for(int iBin = gp->gpProbabilityStat[iRefLayer]->GetNbinsX() + 1; iBin >= 0; iBin--) {//we include overflow bin, because it contains the muons found by the gps with higher pt
        takenMuCnt += gp->gpProbabilityStat[iRefLayer]->GetBinContent(iBin);
        if(takenMuCnt >= targetEff * (double)(allMuCnt) ) {
          float threshold = gp->gpProbabilityStat[iRefLayer]->GetXaxis()->GetBinLowEdge(iBin);
          gp->setThreshold(0, threshold);
          cout<<__FUNCTION__<<":"<<__LINE__<<gp->key()<<" iRefLayer "<<iRefLayer<<" takenMuCnt "<<takenMuCnt
              <<" allMuCnt "<<allMuCnt
              <<" found threshold "<<threshold<<std::endl;
          break;
        }
      }
    }
  }

}


void PatternOptimizer::tuneClassProb(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp) {
  unsigned int iRefHit = omtfCand.getRefHitNumber();
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<omtfCandGp->key()<<" omtfResult\n"<<omtfResult<<std::endl;
  //cout<<__FUNCTION__<<":"<<__LINE__<<" "<<exptCandGp->key()<<" exptResult\n"<<exptResult<<std::endl;

  /**
   * for the event with the muon pt > ptCut
   * for each low pt GP (pt < ptCut) histogramming the difference between the PdfSum() of this GP and the PdfSum of the omtfCandGp
   * if the
   */
  //unsigned int iPdf = omtfConfig->nPdfBins()/2;

  int maxDiffBin = goldenPatterns[goldenPatterns.size()-1]->statisitics[0][0].size();

  //we have to find the gp above the threshold with the highest result
  double maxPdfSum = -1;
  auto& maxResult = goldenPatterns[selectedPatNum]->getResults()[candProcIndx][iRefHit];
  //unsigned int maxResultPat = 0;
  for(unsigned int iPat = selectedPatNum; iPat < goldenPatterns.size(); iPat++) {
    if(goldenPatterns[iPat]->key().thePt == 0)
      continue;

    if(goldenPatterns[iPat]->key().theCharge != goldenPatterns[exptPatNum]->key().theCharge) //exptPatNum should be = selectedPatNum
      continue;

    auto& result = goldenPatterns[iPat]->getResults()[candProcIndx][iRefHit];
    if(result.getFiredLayerCnt() == omtfResult.getFiredLayerCnt() ) {
      if(result.getPdfSum() > maxPdfSum) {
        maxPdfSum = result.getPdfSum();
        maxResult = result;
        //maxResultPat = iPat;
      }
    }
  }

  //if(omtfCandGp->key().thePt >= ptCut) //if commented, taking also the event in which the OMTF cand was below the thresh
  {
    for(unsigned int iPat = 0; iPat < goldenPatterns.size(); iPat++) {
      if(goldenPatterns[iPat]->key().thePt == 0)
        continue;

      if(goldenPatterns[iPat]->key().theCharge != goldenPatterns[exptPatNum]->key().theCharge)
        continue;

      unsigned int iRefLayer = omtfCand.getRefLayer();
      unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];

      if(goldenPatterns[iPat]->key().thePt < ptCut) {
        auto& result = goldenPatterns[iPat]->getResults()[candProcIndx][iRefHit];

        double diff = maxDiffBin -1; //in the bin maxDiffBin -1 will go the events that are below the threshold

        if(maxPdfSum > 0  &&  result.getFiredLayerCnt() == omtfResult.getFiredLayerCnt()) {
          diff = maxResult.getPdfSum() - result.getPdfSum();

          /*if(result.getFiredLayerCnt() == omtfResult.getFiredLayerCnt() &&  diff < 0 ) {//debug
            cout<<__FUNCTION__<<":"<<__LINE__<<" diff = "<<diff<<std::endl;
            cout<<__FUNCTION__<<":"<<__LINE__<<" gp "<<goldenPatterns[iPat]->key()<<"\n"<<result<<std::endl;
            cout<<__FUNCTION__<<":"<<__LINE__<<" maxResultPat "<<goldenPatterns[maxResultPat]->key()<<"\n"<<maxResult<<std::endl;
          }*/

          diff = diff + maxDiffBin/2; //to allow for negative values

          if(diff >= (maxDiffBin-1) ) {
            diff = maxDiffBin -2; //overflows go here
          }
          else if (diff < 0) //watch out - this is already after adding the offset, so it is just to protect against seg fault!
            diff = 0;

          //debug
          /*if(diff < maxDiffBin/2) {
          cout<<__FUNCTION__<<":"<<__LINE__<<" diff < 0 !!!!!! diff = "<<diff<<std::endl;
          cout<<__FUNCTION__<<":"<<__LINE__<<" gp "<<goldenPatterns[iPat]->key()<<"\n"<<result<<std::endl;
          cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCandGp "<<omtfCandGp->key()<<"\n"<<omtfResult<<std::endl;
        }*/
        }
        //cout<<__FUNCTION__<<":"<<__LINE__<<goldenPatterns[iPat]->key()<<" iRefLayer "<<iRefLayer<<" diff = "<<diff<<std::endl;
        goldenPatterns[iPat]->statisitics[refLayerLogicNumber][iRefLayer][diff][ whatExptVal]++; //we use the  iPdf bin to put the diff histogram
      }
    }
  }
}



void PatternOptimizer::tuneClassProb(double targetEffLost) {
  cout<<__FUNCTION__<<":"<<__LINE__<<" targetEff "<<targetEffLost<<std::endl;
  cout<<__FUNCTION__<<": "<<__LINE__<<" Calculating P(C_k) "<<std::endl;
  unsigned int iPdf = omtfConfig->nPdfBins()/2;// <<(omtfConfig->nPdfAddrBits()-1);
  int maxDiffBin = goldenPatterns[goldenPatterns.size()-1]->statisitics[0][0].size();
  
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;

    for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      if(iRefLayer == 0 || iRefLayer == 2) //DT
        targetEffLost = 0.001;
      else if(iRefLayer == 5) //DT
        targetEffLost = 0.001;
      else if(iRefLayer == 1 || iRefLayer == 3 || iRefLayer == 4) //CSC & RE2/3
        targetEffLost = 0.035;
      else if(iRefLayer == 6 || iRefLayer == 7) //bRPC
        targetEffLost = 0.02;

      if(gp->key().thePt < ptCut && gp->key().theCharge == goldenPatterns[selectedPatNum]->key().theCharge) {
        cout<<endl;
        cout<<__FUNCTION__<<":"<<__LINE__<<gp->key()<<" RefLayer "<<iRefLayer<<std::endl;

        double cnt = 0;
        double norm = 0;
        unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];

        for(int iDiff = 0; iDiff < maxDiffBin; iDiff++) {
          cout<<setw(4)<<(iDiff - maxDiffBin/2) <<" ";
        }
        cout<<endl;
        for(int iDiff = 0; iDiff < maxDiffBin; iDiff++) {
          norm += gp->statisitics[refLayerLogicNumber][iRefLayer][iDiff][whatExptVal];
          cout<<setw(4)<<gp->statisitics[refLayerLogicNumber][iRefLayer][iDiff][whatExptVal]<<" ";
        }
        cout<<endl;

        int iDiff = 0;
        for(; iDiff < maxDiffBin -1; iDiff++) { //last but one bin is overflow, last are the events with less fired planes
          cnt += gp->statisitics[refLayerLogicNumber][iRefLayer][iDiff][ whatExptVal];
          if( (cnt / norm) >= targetEffLost) {
            int diff = iDiff - maxDiffBin/2.;
            diff = diff - 1; //if thepdfSum is equal the lower pt GP is chosen, see <class GoldenPatternType> AlgoMuon OMTFSorter<GoldenPatternType>::sortRefHitResults
            cout<<__FUNCTION__<<":"<<__LINE__<<" cnt "<<cnt<<" norm "<<norm<<" current pdfVal "<<gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf]<<" diff "<<diff<<std::endl;

            if(gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf+1] == 0 || diff < 0) {//we are using the iPdf+1 to mark that the value was already updated
              gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf] += diff;

              cout<<__FUNCTION__<<":"<<__LINE__<<" new pdfVal "<<gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf]<<std::endl;
              if(gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf] <= 0) {
                cout<<__FUNCTION__<<":"<<__LINE__<<" pdf <= 0 !!!!!!!!!!!!!!!!!!!!!!"<<endl;;
              }

              gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf+1] = 1;
            }
            break;
          }
        }
        /*if(cnt == 0 && norm > 0)
          diff = maxDiffBin -1;*/
        if(iDiff == maxDiffBin -1) {
          cout<<__FUNCTION__<<":"<<__LINE__<<" diff not found, pdf not updated!!!! cnt "<<cnt<<" norm "<<norm<<" current pdf val "<<gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf]<<" iDiff "<<iDiff<<std::endl;
          //cout<<__FUNCTION__<<":"<<__LINE__<<" assuming diff is 52 !!!!!!!!!!!!!!!!!!!"<<endl;
          //gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf] += 52; ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
      }

    }
  }

/*  for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns[0]->getPdf()[0].size(); ++iRefLayer) {
    for(auto& gp : goldenPatterns) {
      if(gp->key().thePt == 0)
        continue;
      if(gp->key().thePt < ptCut) {

        cout<<__FUNCTION__<<":"<<__LINE__<<gp->key()<<" RefLayer "<<iRefLayer<<std::endl;
      }
    }
  }*/
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
  outfile.mkdir("patternsPdfs/canvases");
  ostringstream ostrName;
  ostringstream ostrTtle;
  vector<TH1F*> classProbHists;
  for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns[0]->getPdf()[0].size(); ++iRefLayer) {
    ostrName.str("");
    ostrName<<"Neg_RefLayer_"<<iRefLayer;
    ostrTtle.str("");
    ostrTtle<<"Neg_RefLayer_"<<iRefLayer;
    classProbHists.push_back(new TH1F(ostrName.str().c_str(), ostrTtle.str().c_str(), goldenPatterns.size(), -0.5, goldenPatterns.size() - 0.5) );

    ostrName.str("");
    ostrName<<"Pos_RefLayer_"<<iRefLayer;
    ostrTtle.str("");
    ostrTtle<<"Pos_RefLayer_"<<iRefLayer;
    classProbHists.push_back(new TH1F(ostrName.str().c_str(), ostrTtle.str().c_str(), goldenPatterns.size(), -0.5, goldenPatterns.size() - 0.5) );
  }
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
    outfile.cd("patternsPdfs");
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
    outfile.cd("patternsPdfs/canvases");
    canvas->Write();
    delete canvas;

    unsigned int iPdf = omtfConfig->nPdfBins()/2;
    for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
      if(gp->key().theCharge == -1) {
        classProbHists[2 * iRefLayer]->Fill(gp->key().theNumber, gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf]);
      }
      else
        classProbHists[2 * iRefLayer +1]->Fill(gp->key().theNumber, gp->pdfAllRef[refLayerLogicNumber][iRefLayer][iPdf]);
    }
  }

  outfile.mkdir("patternsPdfSumStat")->cd();
  for(auto& gp : goldenPatterns) {
    for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      gp->gpProbabilityStat[iRefLayer]->Write();
    }
  }

  outfile.cd();
  for(auto& classProbHist : classProbHists) {
    classProbHist->Write();
  }
  outfile.Close();
}



double vxMuRate(double pt_GeV)
{
  if (pt_GeV == 0)
    return 0.0;
  const double lum = 2.0e34; //defoult is 1.0e34;
  const double dabseta = 1.0;
  const double dpt = 1.0;
  const double afactor = 1.0e-34*lum*dabseta*dpt;
  const double a  = 2*1.3084E6;
  const double mu=-0.725;
  const double sigma=0.4333;
  const double s2=2*sigma*sigma;

  double ptlog10;
  ptlog10 = log10(pt_GeV);
  double ex = (ptlog10-mu)*(ptlog10-mu)/s2;
  double rate = (a * exp(-ex) * afactor);
  //edm::LogError("RPCTrigger")<<ptCode<<" "<<rate;//<<<<<<<<<<<<<<<<<<<<<<<<
  return rate;
 }

double vxIntegMuRate(double pt_GeV, double dpt, double etaFrom, double etaTo) {
  //calkowanie metoda trapezow - nie do konca dobre
  double rate = 0.5 * (vxMuRate(pt_GeV) + vxMuRate(pt_GeV+dpt)) * dpt;

  rate = rate * (etaTo - etaFrom);
  //edm::LogError("RPCTrigger")<<ptCode<<" "<<rate;//<<<<<<<<<<<<<<<<<<<<<<<<
  return rate;
}

void PatternOptimizer::initRateWeights() {
  rateWeights.assign(goldenPatterns.size(), 0);
  patternPtCodes.assign(goldenPatterns.size(), 0);

  int positivePtCode = 1;
  int negativePtCode = 1;
  for(unsigned int patNum = 0; patNum < goldenPatterns.size(); patNum++) {
    if(goldenPatterns[patNum]->key().thePt == 0)
      continue;

    double ptFrom = omtfConfig->getPatternPtRange(patNum).ptFrom;
    double ptTo = omtfConfig->getPatternPtRange(patNum).ptTo;

    rateWeights[patNum] = vxIntegMuRate(ptFrom, ptTo - ptFrom, 0.82, 1.24);

    if(goldenPatterns[patNum]->key().theCharge == 1) {
      patternPtCodes[patNum] = positivePtCode;
      positivePtCode++;
    }
    else if(goldenPatterns[patNum]->key().theCharge == -1) {
      patternPtCodes[patNum] = negativePtCode;
      negativePtCode++;
    }

    std::cout<<goldenPatterns[patNum]->key()<<" ptCode "<<patternPtCodes[patNum]<<" rate "<<rateWeights[patNum]<< std::endl;
  }
}

void PatternOptimizer::modifyPatterns1(double step) {
  cout<<__FUNCTION__<<": "<<__LINE__<<" Correcting P(C_k) "<<std::endl;
  unsigned int iPdf = omtfConfig->nPdfBins()/2;// <<(omtfConfig->nPdfAddrBits()-1);
  double delta = 0;
  for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns[0]->getPdf()[0].size(); ++iRefLayer) {
    unsigned int refLayerLogicNumber = omtfConfig->getRefToLogicNumber()[iRefLayer];
    if(iRefLayer == 0 || iRefLayer == 2) //DT
      step = 1.5;
    else if(iRefLayer == 5) //DT
      step = 1.5;
    else if(iRefLayer == 1) //CSC
      step = 0.33333333;
    else if(iRefLayer == 3) //CSC
          step = 0.5;
    else if(iRefLayer == 5) //RE2/3
      step = 0.5;
    else if(iRefLayer == 6 || iRefLayer == 7) //bRPC
      step = 1.5;

    cout<<__FUNCTION__<<":"<<__LINE__<<" RefLayer "<<iRefLayer<<" step "<<step<<std::endl;
    for(int sign = -1; sign <= 1; sign++) {
      delta = 0;
      for(auto& gp : boost::adaptors::reverse(goldenPatterns) ) {
        if(gp->key().thePt == 0 || gp->key().theCharge != sign)
          continue;
        int newPdfVal = gp->getPdf()[refLayerLogicNumber][iRefLayer][iPdf] + delta;
        gp->setPdfValue(newPdfVal, refLayerLogicNumber, iRefLayer, iPdf);
        delta += step;
      }
    }

  }

}
