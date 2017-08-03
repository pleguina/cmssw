/*
 * PdfGeneratorProcessor.cpp
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PdfGeneratorProcessor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

////////////////////////////////////////////
////////////////////////////////////////////

///////////////////////////////////////////////
///////////////////////////////////////////////
void  PdfGeneratorProcessor::averagePatterns(){
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
  for(unsigned int iGroup = 0; iGroup < mergedPartters.size(); iGroup++) {
    if(mergedPartters[iGroup].size() > 1) {

      GoldenPattern* gp = dynamic_cast<GoldenPattern*>(theGPs.at(mergedPartters[iGroup][0]));
      if(gp == 0)
        throw cms::Exception("PdfGeneratorProcessor::averagePatterns works only for GoldenPatter");

      GoldenPattern::vector3D meanDistPhi = gp->getMeanDistPhi();

      for(unsigned int i = 1; i < mergedPartters[iGroup].size(); i++) {
        GoldenPattern* gp = static_cast<GoldenPattern*>(theGPs.at( mergedPartters[iGroup][i]));
        for(unsigned int iLayer=0; iLayer<myOmtfConfig->nLayers(); ++iLayer) {
          for(unsigned int iRefLayer=0; iRefLayer<myOmtfConfig->nRefLayers(); ++iRefLayer) {
            meanDistPhi[iLayer][iRefLayer][0] += gp->getMeanDistPhi()[iLayer][iRefLayer][0];
          }
        }
      }

      for(unsigned int iLayer=0; iLayer<myOmtfConfig->nLayers(); ++iLayer) {
        for(unsigned int iRefLayer=0; iRefLayer<myOmtfConfig->nRefLayers(); ++iRefLayer) {
          meanDistPhi[iLayer][iRefLayer][0] /= mergedPartters[iGroup].size();
        }
      }

      for(unsigned int i = 0; i < mergedPartters[iGroup].size(); i++) {
        GoldenPattern* gp = static_cast<GoldenPattern*>(theGPs.at( mergedPartters[iGroup][i]));
        shiftGP(gp, meanDistPhi, gp->getMeanDistPhi());
        gp->setMeanDistPhi(meanDistPhi);
      }
    }
  }
}
///////////////////////////////////////////////
///////////////////////////////////////////////
void PdfGeneratorProcessor::shiftGP(GoldenPattern *aGP,
    const GoldenPattern::vector3D & meanDistPhiNew,
    const GoldenPattern::vector3D & meanDistPhiOld){

  ///Shift pdfs by differecne between original menaDistPhi, and
  ///the averaged value
  unsigned int nPdfBins =  exp2(myOmtfConfig->nPdfAddrBits());
  GoldenPattern::vector3D pdfAllRef = aGP->getPdf();

  int indexShift = 0;
  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
    for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
      indexShift = meanDistPhiOld[iLayer][iRefLayer][0] - meanDistPhiNew[iLayer][iRefLayer][0];
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin)
        pdfAllRef[iLayer][iRefLayer][iPdfBin] = 0;
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin){
        if((int)(iPdfBin)+indexShift>=0 && iPdfBin+indexShift<nPdfBins)
          pdfAllRef[iLayer][iRefLayer][iPdfBin+indexShift] = aGP->pdfValue(iLayer, iRefLayer, iPdfBin);
      }
    }
  }
  aGP->setPdf(pdfAllRef);
}
