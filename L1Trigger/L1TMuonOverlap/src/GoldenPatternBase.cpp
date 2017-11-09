/*
 * GoldenPatternBase.cpp
 *
 *  Created on: Oct 3, 2017
 *      Author: kbunkow
 */


#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h"

////////////////////////////////////////////////////
////////////////////////////////////////////////////
GoldenPatternResult::LayerResult GoldenPatternBase::process1Layer1RefLayer(unsigned int iRefLayer,
    unsigned int iLayer,
    const int phiRefHit,
    const OMTFinput::vector1D & layerHits,
    int refLayerPhiB)
{
  //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) std::cout<<__FUNCTION__<<":"<<__LINE__<<key()<<this->getDistPhiBitShift(iLayer, iRefLayer)<<std::endl;
  GoldenPatternResult::LayerResult aResult(0, 0, 0); //0, 0

  int phiMean = this->meanDistPhiValue(iLayer, iRefLayer, refLayerPhiB);
  int phiDistMin = myOmtfConfig->nPdfBins(); //1<<(myOmtfConfig->nPdfAddrBits()); //"infinite" value for the beginning
  ///Select hit closest to the mean of probability
  ///distribution in given layer
  for(auto itHit: layerHits){
    if(itHit >= (int)myOmtfConfig->nPhiBins())
      continue;  //empty itHits are marked with nPhiBins() in OMTFProcessor::restrictInput

    int phiDist = itHit-phiMean-phiRefHit;
    //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) std::cout<<__FUNCTION__<<":"<<__LINE__<<" itHit "<<itHit<<" phiMean "<<phiMean<<" phiRefHit "<<phiRefHit<<" phiDist "<<phiDist<<std::endl;
    phiDist = phiDist >> this->getDistPhiBitShift(iLayer, iRefLayer); //if the shift is done here, it means that the phiMean in the xml should be the same as without shift
    //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) std::cout<<__FUNCTION__<<":"<<__LINE__<<" phiDist "<<phiDist<<std::endl;
    if(abs(phiDist) < abs(phiDistMin))
      phiDistMin = phiDist;
  }

  ///Check if phiDistMin is within pdf range -63 +63
  ///in firmware here the arithmetic "value and sign" is used, therefore the range is -63 +63, and not -64 +63
  if(abs(phiDistMin) > ( (1<<(myOmtfConfig->nPdfAddrBits()-1)) -1) ) {
    //return aResult;
    return GoldenPatternResult::LayerResult(this->pdfValue(iLayer, iRefLayer, 0), false, 0);
    //in the version with thresholds we use the bin 0 to store the pdf value returned when there was no hit.
    //in the version without thresholds, the value in the bin 0 should be 0
  }

  ///Shift phidist, so 0 is at the middle of the range
  phiDistMin += 1<<(myOmtfConfig->nPdfAddrBits()-1);
  //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) std::cout<<__FUNCTION__<<":"<<__LINE__<<" phiDistMin "<<phiDistMin<<std::endl;
  omtfPdfValueType pdfVal = this->pdfValue(iLayer, iRefLayer, phiDistMin);
  if(pdfVal <= 0)
    return GoldenPatternResult::LayerResult(this->pdfValue(iLayer, iRefLayer, 0), false, phiDistMin); //this is needed for the version with threshold
  return GoldenPatternResult::LayerResult(pdfVal, true, phiDistMin);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternBase::finalise() {
  for(auto& result : getResults()) {
    result.finalise();
  }
}
