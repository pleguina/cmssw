#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <cmath>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"

////////////////////////////////////////////
////////////////////////////////////////////
void GoldenPatternResult::configure(const OMTFConfiguration * omtfConfig) {

  myOmtfConfig = omtfConfig;

  reset();

}
////////////////////////////////////////////
////////////////////////////////////////////

void GoldenPatternResult::set(int refLayer_, unsigned int phi, unsigned int eta, unsigned int refHitPhi,
    unsigned int iLayer, unsigned int pdfVal) {
  if( isValid() && this->refLayer != refLayer_) {
    std::cout<<__FUNCTION__<<" "<<__LINE__<<" this->refLayer "<<this->refLayer<<" refLayer_ "<<refLayer_<<std::endl;
  }
  assert( !isValid() || this->refLayer == refLayer_);

  refLayer = refLayer_;
  this->phi = phi;
  this->eta = eta;
  this->refHitPhi = refHitPhi;
  pdfWeights[iLayer] = pdfVal;
  //std::cout<<__FUNCTION__<<" "<<__LINE__<<" iLayer "<<iLayer<<" refLayer "<<refLayer<<std::endl;
  //pdfWeightSum += pdfVal; - this cannot be done here, because the pdfVal for the banding layer must be added only
  //if hit in the corresponding phi layer was accpeted (i.e. its pdfVal > 0. therefore it is done in finalise()
}


/*void GoldenPatternResult::setRefPhiRHits(unsigned int iRefLayer, int iRefPhiRHit) {
  refHitPhi = iRefPhiRHit;
}
////////////////////////////////////////////
////////////////////////////////////////////

void GoldenPatternResult::addResult(unsigned int iRefLayer,
			   unsigned int iLayer,
			   unsigned int val,
			   int iRefPhi, 
			   int iRefEta){

  refLayerResults[iRefLayer].phi = iRefPhi;
  refLayerResults[iRefLayer].eta = iRefEta;
  refLayerResults[iRefLayer].pdfWeights[iLayer] = val;

}*/
////////////////////////////////////////////
////////////////////////////////////////////
void GoldenPatternResult::reset() {
  refLayer = -1;
  pdfWeights.assign(myOmtfConfig->nLayers(), 0);
  phi = 1024;
  eta = 1024;
  pdfWeightSum = 0;
  firedLayerCnt = 0;
  firedLayerBits = 0;
  refHitPhi = 1024;
}

/*void GoldenPatternResult::clear() {
  if(refLayerResults.size() == 0)
    refLayerResults.assign(myOmtfConfig->nRefLayers(), RefLayerResult());
  for (auto& reflayerRes: refLayerResults) {
    reflayerRes.reset();
  }
  results1D.assign(myOmtfConfig->nRefLayers(),0);
  hits1D.assign(myOmtfConfig->nRefLayers(),0);
  results.assign(myOmtfConfig->nLayers(),results1D);
  refPhi1D.assign(myOmtfConfig->nRefLayers(),1024);
  refEta1D.assign(myOmtfConfig->nRefLayers(),1024);
  hitsBits.assign(myOmtfConfig->nRefLayers(),0);  
  refPhiRHit1D.assign(myOmtfConfig->nRefLayers(),1024);
}*/
////////////////////////////////////////////
////////////////////////////////////////////
void GoldenPatternResult::finalise() {
  for(unsigned int iLogicLayer=0; iLogicLayer < pdfWeights.size(); ++iLogicLayer) {
    unsigned int connectedLayer = myOmtfConfig->getLogicToLogic().at(iLogicLayer);

    ///If connected layer (POS or BEND) has not been fired, ignore this layer also
    unsigned int val = pdfWeights[connectedLayer] > 0 ? pdfWeights[iLogicLayer] : 0;
    pdfWeightSum += val;
    firedLayerBits += (val>0)*std::pow(2,iLogicLayer);
    ///Do not count bending layers in hit count
    if(!myOmtfConfig->getBendingLayers().count(iLogicLayer))
      firedLayerCnt += (val>0);
  }
}

////////////////////////////////////////////
////////////////////////////////////////////
/*bool GoldenPatternResults::empty() const{

  unsigned int nHits = 0;
  for(unsigned int iRefLayer=0; iRefLayer<myOmtfConfig->nRefLayers(); ++iRefLayer){
    nHits+=hits1D[iRefLayer];
  }      
  return (nHits==0);
}*/
////////////////////////////////////////////
////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const GoldenPatternResult & gpResult) {
  for(unsigned int iLogicLayer=0; iLogicLayer < gpResult.getPdfWeights().size(); ++iLogicLayer){
    out<<" layer: "<<iLogicLayer<<" res: ";
    out<<std::setw(3)<<gpResult.getPdfWeights()[iLogicLayer]<<" ";
  }
  out<<std::endl;

  out<<"  refLayer: ";
  out << gpResult.getRefLayer()<<"\t";

  out<<" Sum over layers: ";
  out<<gpResult.getPdfWeigtSum()<<"\t";

  out<<" Number of hits: ";
  out << gpResult.getFiredLayerCnt()<<"\t";
  out<<std::endl;


  return out;
}
////////////////////////////////////////////
////////////////////////////////////////////
