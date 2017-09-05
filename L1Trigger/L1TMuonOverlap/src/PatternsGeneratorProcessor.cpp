/*
 * PatternsGeneratorProcessor.cpp
 *
 *  Created on: Aug 3, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PatternsGeneratorProcessor.h"

#include <climits>

PatternsGeneratorProcessor::PatternsGeneratorProcessor(): ProcessorBase<GoldenPattern>() {
  // TODO Auto-generated constructor stub

}

PatternsGeneratorProcessor::~PatternsGeneratorProcessor() {
  // TODO Auto-generated destructor stub
}

void PatternsGeneratorProcessor::generatePatterns() {

}

void PatternsGeneratorProcessor::generateThresholds() {
  for(auto& gp : theGPs) {
    std::vector<unsigned int> thresholds(gp->getPdf()[0].size(), 0);
    for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      std::vector<int> maxPdfValues(gp->getPdf()[0].size());
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
        maxPdfValues[iRefLayer] = *std::max_element(gp->getPdf()[iLayer][iRefLayer].begin(), gp->getPdf()[iLayer][iRefLayer].end() );
        int pdfVal = 0;//maxPdfValues[iRefLayer]
        gp->setPdfValue(0, iLayer, iRefLayer, 0);//for the version with thresholds in the bin 0 of the pdf we  keep the max vale
        //it is used when given layer is not fired
      }

      int  smallestMax  = INT_MAX;
      unsigned int smallestMaxPos = maxPdfValues.size();
      for(unsigned int i = 0; i < maxPdfValues.size(); i++) {
        if(maxPdfValues[i] > 0 && smallestMax > maxPdfValues[i]) {
          smallestMax = maxPdfValues[i];
          smallestMaxPos = i;
        }
      }
      if(smallestMax == INT_MAX) {
        smallestMax = 0;
        smallestMaxPos = maxPdfValues.size();
      }

      thresholds[iRefLayer] = smallestMax;
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
        if(iLayer == smallestMaxPos)
          continue;

        int minPdfVal = INT_MAX;//todo
        for(auto& pdfVal : gp->getPdf()[iLayer][iRefLayer]) {
          if( (pdfVal > 0) && (minPdfVal > pdfVal) )
            minPdfVal = pdfVal;
        }
        if(minPdfVal != INT_MAX)
          thresholds[iRefLayer] += minPdfVal;
        //std::cout<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" minPdfVal "<<minPdfVal<<" thr "<<thresholds[iRefLayer]<<std::endl;
      }
      if(gp->key().thePt <= 10) {
        thresholds[iRefLayer] = 0;
      }
      else if(gp->key().thePt >= 50) {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 1.35;
      }
      else if(gp->key().thePt >= 33) {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 1.35;
      }
      else {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 0.4;
      }
    }
    gp->setThresholds(thresholds);
  }
}
