/*
 * PatternGenerator.cc
 *
 *  Created on: Nov 8, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/PatternGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TFile.h"

PatternGenerator::PatternGenerator(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  PatternOptimizerBase(edmCfg, omtfConfig, gps), eventCntPerGp(gps.size(), 0)
{
  edm::LogImportant("l1tMuBayesEventPrint") << "constructing PatternGenerator " << std::endl;

//TODO uncomment when needed
/*  //adding new patterns!!!!!!!!!!!!
  edm::LogImportant("PatternGeneratorTT") << "PatternGeneratorTT: adding new patterns and modifying existing!!!!!" << std::endl;

  gps[46]->setKeyPt(45);
  gps[47]->setKeyPt(49);
  gps[52]->setKeyPt(53);

  gps[42]->setKeyPt(45);
  gps[43]->setKeyPt(49);
  gps[48]->setKeyPt(53);

  todo 0 set also group and indexInGroup, in any case not easy
  Easier is to add the new patterns in the input xml, just set the iPt3 or iPt4 accordingly

  auto pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 0, -1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 0, -1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 8, -1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 7, -1, 0), omtfConfig)); pos = gps.begin();

  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 0, 1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 0, 1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 8, 1, 0), omtfConfig)); pos = gps.begin();
  gps.insert(pos, make_shared<GoldenPatternWithStat>(Key(0, 7, 1, 0), omtfConfig)); pos = gps.begin(); */

  goldenPatterns = gps;

  //reseting the golden patterns
  unsigned int i = 0;
  for(auto& gp : goldenPatterns) {
      gp->setKeyNumber(i++); //needed  if patterns were added

      if(gp->key().thePt == 0)
        continue;

      gp->reset();

      int statBinsCnt = 1024; //gp->getPdf()[0][0].size() * 8; //TODO should be big enough to comprise the pdf tails
      gp->iniStatisitics(statBinsCnt, 1); //TODO
  }

  //GoldenPatternResult::setFinalizeFunction(3); TODO why it was this one????
 // edm::LogImportant("l1tMuBayesEventPrint") << "reseting golden pattern !!!!!" << std::endl;

  //setting all pdf to 1, this will cause that the when the OmtfProcessor process the input, the result will be based only on the number of fired layers,
  //and then the omtfCand will come from the processor that has the biggest number of fired layers
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        //unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
        //if(refLayerLogicNum == iLayer)
        {
          for(unsigned int iBin = 0; iBin < gp->getPdf()[iLayer][iRefLayer].size(); iBin++) {
            gp->pdfAllRef[iLayer][iRefLayer][iBin] = 1;
          }
        }
      }
    }
  }

/*  ptDeltaPhiHists.resize(2);
  for(unsigned int iCharge = 0; iCharge <= 1; iCharge++) {
    for(unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) { //for the moment filing only ref layer, remove whe
      if(iLayer == 0 || iLayer == 2 || iLayer == 4 || iLayer == 6 || iLayer == 7 || iLayer == 10 || iLayer == 11  || iLayer == 16 || //refLayars
         iLayer == 1 || iLayer == 3 || iLayer == 5  ) //banding layers
      {
        ostringstream name;
        name<<"ptDeltaPhiHist_ch_"<<iCharge<<"_Layer_"<<iLayer;
        int phiFrom = -10;
        int phiTo   = 300; //TODO
        int phiBins = phiTo - phiFrom;

        if(iCharge == 1) {
          phiFrom = -300; //TODO
          phiTo = 10;
        }

        TH2I* ptDeltaPhiHist = new TH2I(name.str().c_str(), name.str().c_str(), 400, 0, 200, phiBins, phiFrom -0.5, phiTo -0.5);
        //cout<<"BinLowEdge "<<ptDeltaPhiHist->GetYaxis()->GetBinLowEdge(100)<<" BinUpEdge "<<ptDeltaPhiHist->GetYaxis()->GetBinUpEdge(100);
        ptDeltaPhiHists[iCharge].push_back(ptDeltaPhiHist);
      }
      else
        ptDeltaPhiHists[iCharge].push_back(nullptr);
    }
  }*/
}

PatternGenerator::~PatternGenerator() {

}

void PatternGenerator::updateStat() {
  //cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCand "<<*omtfCand<<std::endl;;
  AlgoMuon* algoMuon = omtfCand.get();
  if(!algoMuon) {
    edm::LogImportant("l1tMuBayesEventPrint")<<":"<<__LINE__<<" algoMuon is null"<<std::endl;
    throw runtime_error("algoMuon is null");
  }

  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;

  unsigned int exptPatNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(exptPatNum).get(); // expected pattern

  eventCntPerGp[exptPatNum]++;

  //edm::LogImportant("l1tMuBayesEventPrint")<<"\n" <<__FUNCTION__<<": "<<__LINE__<<" exptCandGp "<<exptCandGp->key()<<" candProcIndx "<<candProcIndx<<" ptSim "<<ptSim<<" chargeSim "<<chargeSim<<std::endl;

  unsigned int iCharge = omtfCand->getCharge();
  if(iCharge != 1)
    iCharge = 0;

  int pdfMiddle = 1<<(omtfConfig->nPdfAddrBits()-1);

 //iRefHit is the index of the hit
  for(unsigned int iRefHit = 0; iRefHit < exptCandGp->getResults()[candProcIndx].size(); ++iRefHit) {
    auto& gpResult  = exptCandGp->getResults()[candProcIndx][iRefHit];
    //unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefHit];

    unsigned int refLayer = gpResult.getRefLayer();

    //cout<<gpResult;
    if(gpResult.getFiredLayerCnt() >= 3 )
    {
      //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic"<<std::endl;
      for(unsigned int iLayer = 0;  iLayer < gpResult.getStubResults().size(); iLayer++) {
        //updating statistic for the gp which should have fired
        if(gpResult.getStubResults()[iLayer].getMuonStub() ) {//the result is not empty
          int phiDist = gpResult.getStubResults()[iLayer].getPdfBin();
          phiDist += exptCandGp->meanDistPhiValue(iLayer, refLayer) - pdfMiddle; //removing the shift applied in the GoldenPatternBase::process1Layer1RefLayer

          /*
          if(ptDeltaPhiHists[iCharge][iLayer] != nullptr &&
              (iLayer == refLayerLogicNum || omtfConfig->getLogicToLogic().at(iLayer) == (int)refLayerLogicNum) )
            ptDeltaPhiHists[iCharge][iLayer]->Fill(ttAlgoMuon->getPt(), phiDist); //TODO correct
           */

          phiDist += exptCandGp->getStatistics()[iLayer][refLayer].size()/2;

          //edm::LogImportant("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" refLayer "<<refLayer<<" iLayer "<<iLayer<<" phiDist "<<phiDist<<" getPdfBin "<<gpResult.getStubResults()[iLayer].getPdfBin()<<std::endl;
          if( phiDist >= 0 && phiDist < (int)(exptCandGp->getStatistics()[iLayer][refLayer].size()) ) {
            //updating statistic for the gp which found the candidate
            //edm::LogImportant("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" updating statistic "<<std::endl;
            exptCandGp->updateStat(iLayer, refLayer, phiDist, 0, 1);
          }
        }
      }
    }
  }

}

void PatternGenerator::observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  if(simMuon == 0 || omtfCand->getGoldenPatern() == 0)//no sim muon or empty candidate
    return;

  PatternOptimizerBase::observeEventEnd(iEvent, finalCandidates);

  updateStat();
}

void PatternGenerator::endJob() {

  upadatePdfs();

  PatternOptimizerBase::endJob();
}

void PatternGenerator::upadatePdfs() {
  //TODO setting the DistPhiBitShift i.e. grouping of the pdfBins
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        if(gp->getDistPhiBitShift(iLayer, iRefLayer) ) {
          throw runtime_error(string(__FUNCTION__) + ":" + to_string(__LINE__) + "gp->getDistPhiBitShift(iLayer, iRefLayer) != 0 -  cannot change DistPhiBitShift then!!!!");
        }

        if( (gp->key().thePt <= 10) && (iLayer < 6) ) {
          gp->setDistPhiBitShift(1, iLayer, iRefLayer);
        }
        else
          gp->setDistPhiBitShift(0, iLayer, iRefLayer);

        //watch out: the shift in a given layer must be the same for patterns in one group
        //todo  make the setting on shift on the group base
        /*if( (gp->key().thePt <= 10) && (iLayer == 3 || iLayer == 5 ) && (iRefLayer == 0 || iRefLayer == 2 || iRefLayer == 6 || iRefLayer == 7)) {
          gp->setDistPhiBitShift(3, iLayer, iRefLayer);
        }
        else if( (gp->key().thePt <= 10) && ( iLayer == 1 || iLayer == 3 || iLayer == 5 ) ) {
          gp->setDistPhiBitShift(2, iLayer, iRefLayer);
        }
        else if( ( (gp->key().thePt <= 10) && (iLayer == 7 ||iLayer == 8 || iLayer == 17 ) ) ) {
          gp->setDistPhiBitShift(1, iLayer, iRefLayer);
        }
        else if( (gp->key().thePt <= 10) && (iLayer == 10 || iLayer == 11 || iLayer == 12 || iLayer == 13) && (iRefLayer == 1)) {
          gp->setDistPhiBitShift(1, iLayer, iRefLayer);
        }*/
      }
    }
  }


  //Calculating meanDistPhi
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;

    edm::LogImportant("l1tMuBayesEventPrint") << "PatternGenerator::upadatePdfs() "<<gp->key()<<" eventCnt "<<eventCntPerGp[gp->key().number()] << std::endl;
    int minHitCnt = 0.001 * eventCntPerGp[gp->key().number()];// //TODO tune threshold <<<<<<<<<<<<<<<<<<
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        //unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
        //if(refLayerLogicNum == iLayer)
        {
          //calculate meanDistPhi
          double meanDistPhi = 0;
          double count = 0;
          for(unsigned int iBin = 0; iBin < gp->getStatistics()[iLayer][iRefLayer].size(); iBin++) {
            meanDistPhi +=  iBin * gp->getStatistics()[iLayer][iRefLayer][iBin][0];
            count       +=         gp->getStatistics()[iLayer][iRefLayer][iBin][0];
          }

          if(count != 0) {
            meanDistPhi /= count;

            meanDistPhi -= (gp->getStatistics()[iLayer][iRefLayer].size() / 2);

            if(count < minHitCnt)
              meanDistPhi = 0;
            else
              edm::LogImportant("l1tMuBayesEventPrint") <<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" count "<<count<<" meanDistPhi "<<meanDistPhi<<endl;
          }
          gp->setMeanDistPhiValue(round(meanDistPhi), iLayer, iRefLayer);
        }
      }
    }
  }

  OMTFConfiguration::vector2D patternGroups = omtfConfig->getPatternGroups(goldenPatterns);
  edm::LogImportant("l1tMuBayesEventPrint") <<"patternGroups:"<<std::endl;
  for(unsigned int iGroup = 0; iGroup < patternGroups.size(); iGroup++) {
    edm::LogImportant("l1tMuBayesEventPrint") <<"patternGroup "<<std::setw(2)<<iGroup<<" ";
    for(unsigned int i = 0; i < patternGroups[iGroup].size(); i++) {
      edm::LogImportant("l1tMuBayesEventPrint")<<i<<" patNum "<<patternGroups[iGroup][i]<<" ";
    }
    edm::LogImportant("l1tMuBayesEventPrint")<<std::endl;
  }

//averaging the meanDistPhi for the gp belonging to the same group
  for(unsigned int iLayer = 0; iLayer < goldenPatterns.at(0)->getPdf().size(); ++iLayer) {
    for(unsigned int iRefLayer = 0; iRefLayer < goldenPatterns.at(0)->getPdf()[iLayer].size(); ++iRefLayer) {
      //unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
      //if(refLayerLogicNum == iLayer)
      {
        for(unsigned int iGroup = 0; iGroup < patternGroups.size(); iGroup++) {
          double meanDistPhi = 0;
          int mergedCnt = 0;
          for(unsigned int i = 0; i < patternGroups[iGroup].size(); i++) {
            auto gp = goldenPatterns.at(patternGroups[iGroup][i]).get();
            meanDistPhi += gp->meanDistPhiValue(iLayer, iRefLayer);
            if(gp->meanDistPhiValue(iLayer, iRefLayer) != 0)
              mergedCnt++;
          }

          if(mergedCnt) {
            meanDistPhi /= mergedCnt; //because for some gps the statistics can be too low, and then the meanDistPhiValue is 0, so it should not contribute
            for(unsigned int i = 0; i < patternGroups[iGroup].size(); i++) {
              auto gp = goldenPatterns.at(patternGroups[iGroup][i]).get();
              gp->setMeanDistPhiValue(round(meanDistPhi), iLayer, iRefLayer);
              edm::LogImportant("l1tMuBayesEventPrint") <<__FUNCTION__<<": "<<__LINE__<<" iGroup "<<iGroup<<" numInGroup "<<i<<" "<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" meanDistPhi after averaging "<<meanDistPhi<<endl;
            }
          }
        }
      }
    }
  }


  //calculating the pdfs
  for(auto& gp : goldenPatterns) {
    if(gp->key().thePt == 0)
      continue;
    int minHitCnt = 0.001 * eventCntPerGp[gp->key().number()];// //TODO tune threshold <<<<<<<<<<<<<<<<<<

    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        //unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
        //if(refLayerLogicNum == iLayer)
        {
          double norm = 0;
          for(unsigned int iBin = 0; iBin < gp->getStatistics()[iLayer][iRefLayer].size(); iBin++) {
            norm += gp->getStatistics()[iLayer][iRefLayer][iBin][0];
          }

          int pdfMiddle = gp->getPdf()[iLayer][iRefLayer].size() / 2;
          int statBinGroupSize = 1<<gp->getDistPhiBitShift(iLayer, iRefLayer);
          for(unsigned int iBinPdf = 0; iBinPdf < gp->getPdf()[iLayer][iRefLayer].size(); iBinPdf++) {
            double pdfVal = 0;
            int groupedBins = 0;
            for(int i = 0; i < statBinGroupSize; i++) {
              int iBinStat = statBinGroupSize * ((int)(iBinPdf) - pdfMiddle) + i + gp->meanDistPhiValue(iLayer, iRefLayer);

              iBinStat += (gp->getStatistics()[iLayer][iRefLayer].size()/2);

              if(iBinStat >= 0 && iBinStat < (int)gp->getStatistics()[iLayer][iRefLayer].size() ) {
                pdfVal += gp->getStatistics()[iLayer][iRefLayer][iBinStat][0];
                groupedBins++;
                //cout<<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<" iLayer "<<iLayer<<" iBinStat "<<iBinStat<<" iBinPdf "<<iBinPdf<<" statVal "<<gp->getStatistics()[iLayer][iRefLayer][iBinStat][0]<<endl;
              }
            }
            if(norm > minHitCnt) {
              pdfVal /= norm;
            }
            else
              pdfVal = 0;

            double minPdfValFactor = 1.;
            const double minPlog =  log(omtfConfig->minPdfVal() * minPdfValFactor);
            const double pdfMaxVal = omtfConfig->pdfMaxValue();

            int digitisedVal = 0;
            if(pdfVal >= omtfConfig->minPdfVal() * minPdfValFactor) {
              digitisedVal = rint(pdfMaxVal - log(pdfVal) / minPlog * pdfMaxVal);
            }

            gp->setPdfValue(digitisedVal, iLayer, iRefLayer, iBinPdf );
            //cout<<__FUNCTION__<<": "<<__LINE__<<" "<<gp->key()<<" iLayer "<<iLayer<<" iBinPdf "<<iBinPdf<<" pdfVal "<<pdfVal<<" digitisedVal "<<digitisedVal<<endl;
          }

        }
      }
    }
  }
}

void PatternGenerator::saveHists(TFile& outfile) {
  outfile.mkdir("ptDeltaPhiHists" )->cd();
/*  for(unsigned int iCharge = 0; iCharge <= 1; iCharge++) {
    for(unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) { //for the moment filing only ref layer, remove whe
      if(ptDeltaPhiHists[iCharge][iLayer]) {
        ptDeltaPhiHists[iCharge][iLayer]->Write();
      }
    }
  }*/

  //TODO rather store the gp->getStatistics(), so the pdf calcualtion can be done without runnign on the data every time
}

