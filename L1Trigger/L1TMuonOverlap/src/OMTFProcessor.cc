#include <iostream>
#include <algorithm>
#include <strstream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlap/interface/GhostBuster.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
///////////////////////////////////////////////
///////////////////////////////////////////////
OMTFProcessor::~OMTFProcessor(){

  for(auto it: theGPs) delete it;

}
///////////////////////////////////////////////
///////////////////////////////////////////////
void OMTFProcessor::resetConfiguration(){

  //myResults.clear();
  for(auto it: theGPs) delete it;
  theGPs.clear();
}
///////////////////////////////////////////////
///////////////////////////////////////////////
bool OMTFProcessor::configure(const OMTFConfiguration * omtfConfig,
    const L1TMuonOverlapParams * omtfPatterns){

  resetConfiguration();

  myOmtfConfig = omtfConfig;

  setSorter(new OMTFSorter()); //initialize with the default sorter
  setGhostBuster(new GhostBuster()); //initialize with the default sorter

  //myResults.assign(myOmtfConfig->nTestRefHits(),OMTFProcessor::resultsMap());

  const l1t::LUT* chargeLUT =  omtfPatterns->chargeLUT();
  const l1t::LUT* etaLUT =  omtfPatterns->etaLUT();
  const l1t::LUT* ptLUT =  omtfPatterns->ptLUT();
  const l1t::LUT* pdfLUT =  omtfPatterns->pdfLUT();
  const l1t::LUT* meanDistPhiLUT =  omtfPatterns->meanDistPhiLUT();

  unsigned int nGPs = myOmtfConfig->nGoldenPatterns();
  unsigned int address = 0;
  unsigned int iEta, iPt;
  int iCharge;
  for(unsigned int iGP=0;iGP<nGPs;++iGP){
    address = iGP;
    iEta = etaLUT->data(address);
    iCharge = chargeLUT->data(address)==0? -1:1;
    iPt = ptLUT->data(address);

    GoldenPattern::vector2D meanDistPhi2D(myOmtfConfig->nLayers());
    GoldenPattern::vector1D pdf1D(exp2(myOmtfConfig->nPdfAddrBits()));
    GoldenPattern::vector3D pdf3D(myOmtfConfig->nLayers());
    GoldenPattern::vector2D pdf2D(myOmtfConfig->nRefLayers());
    ///Mean dist phi data
    for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
      GoldenPattern::vector1D meanDistPhi1D(myOmtfConfig->nRefLayers());
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        address = iRefLayer + iLayer*myOmtfConfig->nRefLayers() + iGP*(myOmtfConfig->nRefLayers()*myOmtfConfig->nLayers());
        meanDistPhi1D[iRefLayer] = meanDistPhiLUT->data(address) - (1<<(meanDistPhiLUT->nrBitsData() -1));
      }
      meanDistPhi2D[iLayer] = meanDistPhi1D;    
      ///Pdf data
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        pdf1D.assign(1<<myOmtfConfig->nPdfAddrBits(),0);
        for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf){
          address = iPdf + iRefLayer*(1<<myOmtfConfig->nPdfAddrBits()) +
              iLayer*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits()) +
              iGP*myOmtfConfig->nLayers()*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits());
          pdf1D[iPdf] = pdfLUT->data(address);
        }
        pdf2D[iRefLayer] = pdf1D;
      }
      pdf3D[iLayer] = pdf2D;
    }
    Key aKey(iEta,iPt,iCharge,iGP);

    GoldenPattern *aGP = new GoldenPattern(aKey, myOmtfConfig);
    aGP->setMeanDistPhi(meanDistPhi2D);
    aGP->setPdf(pdf3D);
    addGP(aGP);    
  }
  return true;
}
///////////////////////////////////////////////
///////////////////////////////////////////////
bool OMTFProcessor::addGP(IGoldenPattern *aGP){
/*  auto gpIt = std::find(std::begin(theGPs), std::end(theGPs), aGP->key());
  if(gpIt !=theGPs.end()){
    throw cms::Exception("Corrupted Golden Patterns data")
    <<"OMTFProcessor::addGP(...) "
    <<" Reading two Golden Patterns with the same key: "
    <<aGP->key()<<std::endl;
  }*/

  //FIXME - should we also check the pt, charge and eta???
  if(theGPs.size() != aGP->key().theNumber) {
    throw cms::Exception("Corrupted Golden Patterns data")
    <<"OMTFProcessor::addGP(...) "
    <<" theGPs.size() != aGP->key().theNumber: "
    <<aGP->key()<<std::endl;
  }
  //else theGPs[aGP->key()] = new GoldenPattern(*aGP);
  else theGPs.push_back(aGP);

  for(auto & itResult: aGP->getResults()){
    itResult.configure(myOmtfConfig);
  }

  return true;
}
///////////////////////////////////////////////
///////////////////////////////////////////////
void  OMTFProcessor::averagePatterns(int charge){
/* FIXME
  Key aKey(0, 9, charge);

  while(theGPs.find(aKey)!=theGPs.end()){

    GoldenPattern *aGP1 = (GoldenPattern*)(theGPs.find(aKey)->second);
    GoldenPattern *aGP2 = aGP1;
    GoldenPattern *aGP3 = aGP1;
    GoldenPattern *aGP4 = aGP1;

    ++aKey.thePtCode;
    while(theGPs.find(aKey)==theGPs.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
    if(aKey.thePtCode<=401 && theGPs.find(aKey)!=theGPs.end()) aGP2 =  (GoldenPattern*)(theGPs.find(aKey)->second);

    if(aKey.thePtCode>71){
      ++aKey.thePtCode;
      while(theGPs.find(aKey)==theGPs.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
      if(aKey.thePtCode<=401 && theGPs.find(aKey)!=theGPs.end()) aGP3 =  (GoldenPattern*)(theGPs.find(aKey)->second);

      ++aKey.thePtCode;
      while(theGPs.find(aKey)==theGPs.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
      if(aKey.thePtCode<=401 && theGPs.find(aKey)!=theGPs.end()) aGP4 =  (GoldenPattern*)(theGPs.find(aKey)->second);
    }
    else{
      aGP3 = aGP1;
      aGP4 = aGP2;
    }    
    //HACK. Have to clean this up.
    ///Previously pt codes were going by steps of 1, now this is not the case
    ++aKey.thePtCode;
    while(theGPs.find(aKey)==theGPs.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
    ///////////////////////////////
    
    
    GoldenPattern::vector2D meanDistPhi  = aGP1->getMeanDistPhi();

    GoldenPattern::vector2D meanDistPhi1  = aGP1->getMeanDistPhi();
    GoldenPattern::vector2D meanDistPhi2  = aGP2->getMeanDistPhi();
    GoldenPattern::vector2D meanDistPhi3  = aGP3->getMeanDistPhi();
    GoldenPattern::vector2D meanDistPhi4  = aGP4->getMeanDistPhi();
   
    for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        meanDistPhi[iLayer][iRefLayer]+=meanDistPhi2[iLayer][iRefLayer];
        meanDistPhi[iLayer][iRefLayer]+=meanDistPhi3[iLayer][iRefLayer];
        meanDistPhi[iLayer][iRefLayer]+=meanDistPhi4[iLayer][iRefLayer];
        meanDistPhi[iLayer][iRefLayer]/=4;
      }
    }
    
    aGP1->setMeanDistPhi(meanDistPhi);
    aGP2->setMeanDistPhi(meanDistPhi);

    shiftGP(aGP1,meanDistPhi, meanDistPhi1);
    shiftGP(aGP2,meanDistPhi, meanDistPhi2);   
    if(aGP3!=aGP1 && aGP4!=aGP2){
      aGP3->setMeanDistPhi(meanDistPhi);
      aGP4->setMeanDistPhi(meanDistPhi);
      shiftGP(aGP3,meanDistPhi, meanDistPhi3);   
      shiftGP(aGP4,meanDistPhi, meanDistPhi4);   
    }
  }*/
}
///////////////////////////////////////////////
///////////////////////////////////////////////
void OMTFProcessor::shiftGP(GoldenPattern *aGP,
    const GoldenPattern::vector2D & meanDistPhiNew,
    const GoldenPattern::vector2D & meanDistPhiOld){

  ///Shift pdfs by differecne between original menaDistPhi, and
  ///the averaged value
  unsigned int nPdfBins =  exp2(myOmtfConfig->nPdfAddrBits());
  GoldenPattern::vector3D pdfAllRef = aGP->getPdf();

  int indexShift = 0;
  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
    for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
      indexShift = meanDistPhiOld[iLayer][iRefLayer] - meanDistPhiNew[iLayer][iRefLayer];
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin) pdfAllRef[iLayer][iRefLayer][iPdfBin] = 0;
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin){
        if((int)(iPdfBin)+indexShift>=0 && iPdfBin+indexShift<nPdfBins)
          pdfAllRef[iLayer][iRefLayer][iPdfBin+indexShift] = aGP->pdfValue(iLayer, iRefLayer, iPdfBin);
      }
    }
  }
  aGP->setPdf(pdfAllRef);
}
///////////////////////////////////////////////
///////////////////////////////////////////////
const std::vector<IGoldenPattern*> & OMTFProcessor::getPatterns() const{ return theGPs; }
///////////////////////////////////////////////
///////////////////////////////////////////////
//const std::vector<OMTFProcessor::resultsMap> &
const void OMTFProcessor::processInput(unsigned int iProcessor,
    const OMTFinput & aInput){
/*  for(auto & itRegion: myResults)
  	for(auto & itKey: itRegion) 
  		itKey.second.reset();*/

  for(auto& itGP: theGPs) {
    for(auto& result : itGP->getResults()) {
      result.reset();
    }
  }

  //////////////////////////////////////
  //////////////////////////////////////  
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if(refHitsBits.none())
    return; // myResults;

  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer) {
    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);

    if(!layerHits.size()) continue;
    ///Number of reference hits to be checked. 
    unsigned int nTestedRefHits = myOmtfConfig->nTestRefHits();
    for(unsigned int iRefHit = 0; iRefHit < myOmtfConfig->nRefHits(); ++iRefHit) { //loop over all possible refHits, i.e. 128
      if(!refHitsBits[iRefHit]) continue;
      if(nTestedRefHits-- == 0) break;

      const RefHitDef& aRefHitDef = myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int phiRef = aInput.getLayerData(myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer])[aRefHitDef.iInput];
      int etaRef = aInput.getLayerData(myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer],true)[aRefHitDef.iInput];
      unsigned int iRegion = aRefHitDef.iRegion;

      if(myOmtfConfig->getBendingLayers().count(iLayer)) //this iLayer is a banding layer
        phiRef = 0;  //then in the delta_phi one obtains simply the iLayer phi

      const OMTFinput::vector1D restrictedLayerHits = restrictInput(iProcessor, iRegion, iLayer,layerHits);
      //std::cout<<"iLayer "<<iLayer<<" refHitNum "<<myOmtfConfig->nTestRefHits()-nTestedRefHits-1<<" iRefHit "<<iRefHit;
      //std::cout<<" nTestedRefHits "<<nTestedRefHits<<" aRefHitDef "<<aRefHitDef<<std::endl;

      for(auto& itGP: theGPs) {
        GoldenPattern::layerResult aLayerResult = itGP->process1Layer1RefLayer(aRefHitDef.iRefLayer,iLayer,
            phiRef,
            restrictedLayerHits);
        int phiRefSt2 = itGP->propagateRefPhi(phiRef, etaRef, aRefHitDef.iRefLayer);
/*        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].setRefPhiRHits(aRefHitDef.iRefLayer, phiRef);
        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].addResult(aRefHitDef.iRefLayer, iLayer,
            aLayerResult.first,
            phiRefSt2, etaRef);*/

        //myResults.at(myOmtfConfig->nTestRefHits()-nTestedRefHits-1).at(itGP->key()).set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef, iLayer, aLayerResult.first);
        itGP->getResults().at(myOmtfConfig->nTestRefHits()-nTestedRefHits-1).set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef, iLayer, aLayerResult.first);
      }
    }
  }
  //////////////////////////////////////
  ////////////////////////////////////// 
  {
    /*unsigned int iRefHitNum  = 0;
    for(auto & itRefHit: myResults) {
      for(auto & itKey: itRefHit) {
        itKey.second.finalise();
        if(itKey.second.isValid())
          std::cout<<"iRefHitNum "<<iRefHitNum<<" "<<itKey.first<<"\t"<<itKey.second<<std::endl;
      }
      iRefHitNum++;
    }*/

    for(auto& itGP: theGPs) {
      itGP->finalise();
    }
  }

  std::ostringstream myStr;
  myStr<<"iProcessor: "<<iProcessor<<std::endl;
  myStr<<"Input: ------------"<<std::endl;
  myStr<<aInput<<std::endl; 
  edm::LogInfo("OMTF processor")<<myStr.str();


  return; //myResults;
}   
////////////////////////////////////////////
////////////////////////////////////////////
OMTFinput::vector1D OMTFProcessor::restrictInput(unsigned int iProcessor,
    unsigned int iRegion,
    unsigned int iLayer,
    const OMTFinput::vector1D & layerHits){

  OMTFinput::vector1D myHits = layerHits;

  unsigned int iStart = myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].first;
  unsigned int iEnd = iStart + myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].second -1;

  for(unsigned int iInput=0;iInput<14;++iInput){    
    if(iInput<iStart || iInput>iEnd) myHits[iInput] = myOmtfConfig->nPhiBins();
  }  
  return myHits;
}
////////////////////////////////////////////
////////////////////////////////////////////
void OMTFProcessor::fillCounts(unsigned int iProcessor,
    const OMTFinput & aInput,
    const SimTrack* aSimMuon) {

  int theCharge = (abs(aSimMuon->type()) == 13) ? aSimMuon->type()/-13 : 0; 
  unsigned int  iPt =  RPCConst::iptFromPt(aSimMuon->momentum().pt());
  ///Stupid conersion. Have to go through PAC pt scale, as we later
  ///shift resulting pt code by +1
  iPt+=1;
  if(iPt>31) iPt=200*2+1;
  else iPt = RPCConst::ptFromIpt(iPt)*2.0+1;//MicroGMT has 0.5 GeV step size, with lower bin edge  (uGMT_pt_code - 1)*step_size
  //////

  //////////////////////////////////////  
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if(refHitsBits.none())
    return;

  std::ostringstream myStr;
  myStr<<"iProcessor: "<<iProcessor<<std::endl;
  myStr<<"Input: ------------"<<std::endl;
  myStr<<aInput<<std::endl;
  edm::LogInfo("OMTF processor")<<myStr.str();

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" muon iPt "<<iPt<<" theCharge "<<theCharge<<std::endl;
  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){

    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);
    if(!layerHits.size())
      continue;

    ///Number of reference hits to be checked. 
    for(unsigned int iRefHit=0;iRefHit<myOmtfConfig->nRefHits();++iRefHit){

      if(!refHitsBits[iRefHit]) continue;
      const RefHitDef & aRefHitDef = myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      int phiRef = aInput.getLayerData(refLayerLogicNumber)[aRefHitDef.iInput];

      int refLayerPhiB = 0;
      if(refLayerLogicNumber < 6) {//is DT layer TODO - check
        refLayerPhiB = aInput.getLayerData(refLayerLogicNumber+1)[aRefHitDef.iInput];
        if(myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) == 0) {
          throw cms::Exception("not good: the layer is not bending layer");
        }
      }


      unsigned int iRegion = aRefHitDef.iRegion;
      if(myOmtfConfig->getBendingLayers().count(iLayer))
        phiRef = 0;
      const OMTFinput::vector1D restrictedLayerHits = restrictInput(iProcessor, iRegion, iLayer,layerHits);
      for(auto& itGP: theGPs){
/*        if(itGP->key().theCharge != theCharge)  //TODO it was in the orginal code, i commented it out, check why
          continue;*/
        if(itGP->key().thePtCode != iPt)
          continue;
        itGP->addCount(aRefHitDef.iRefLayer, iLayer, phiRef, restrictedLayerHits, refLayerPhiB);
        //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iLayer "<<iLayer<<" refLayerPhiB "<<refLayerPhiB<<std::endl;
      }
    }
  }
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
bool OMTFProcessor::checkHitPatternValidity(unsigned int hits){

  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<unsigned int> badPatterns = {99840, 34304, 3075, 36928, 12300, 98816, 98944, 33408, 66688, 66176, 7171, 20528, 33856, 35840, 4156, 34880};

  for(auto aHitPattern: badPatterns){
    if(hits==aHitPattern) return false;
  }

  return true;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

std::vector<AlgoMuon> OMTFProcessor::sortResults(int charge) {
  std::vector<AlgoMuon> algoCandidates;
  sorter->sortResults(getPatterns(), algoCandidates, charge);
  return algoCandidates;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
std::vector<l1t::RegionalMuonCand> OMTFProcessor::getFinalcandidates(unsigned int iProcessor, l1t::tftype mtfType, const std::vector<AlgoMuon> & algoCands)
{

  std::vector<l1t::RegionalMuonCand> result;

  for(auto myCand: algoCands){
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(myCand.getPt());
    candidate.setHwEta(myCand.getEta());

    int phiValue = myCand.getPhi();
    if(phiValue>= int(myOmtfConfig->nPhiBins()) )
      phiValue -= myOmtfConfig->nPhiBins();
    ///conversion factor from OMTF to uGMT scale: 5400/576
//    phiValue/=9.375;
    phiValue *= (437./pow(2,12));    // ie. use as in hw: 9.3729977
    candidate.setHwPhi(phiValue);

    candidate.setHwSign(myCand.getCharge()<0 ? 1:0  );
    candidate.setHwSignValid(1);

    unsigned int quality = checkHitPatternValidity(myCand.getFiredLayerBits()) ? 0 | (1 << 2) | (1 << 3)
                                                                     : 0 | (1 << 2);
    if (    abs(myCand.getEta()) == 115
        && (    static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001100000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001010000000").to_ulong()
           )
       ) quality =4;

//  if (abs(myCand.getEta()) == 121) quality = 4;
    if (abs(myCand.getEta()) == 121) quality = 0; // changed on request from HI

    candidate.setHwQual (quality);

    std::map<int, int> trackAddr;
    trackAddr[0] = myCand.getFiredLayerBits();
    trackAddr[1] = myCand.getRefLayer();
    trackAddr[2] = myCand.getDisc();
    candidate.setTrackAddress(trackAddr);
    candidate.setTFIdentifiers(iProcessor,mtfType);
    if (candidate.hwPt() >= 0)  result.push_back(candidate);
  }
  return result;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

