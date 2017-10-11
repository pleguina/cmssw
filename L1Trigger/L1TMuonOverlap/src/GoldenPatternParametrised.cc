#include <iostream>
#include <iomanip>
#include <cmath>

#include <TF1.h>
#include <TLinearFitter.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternParametrised.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"


////////////////////////////////////////////////////
////////////////////////////////////////////////////
GoldenPatternParametrised::GoldenPatternParametrised(const Key& aKey, const OMTFConfiguration * omtfConfig) : GoldenPattern(aKey, omtfConfig) {
  //reset(); called by GoldenPattern()
}

GoldenPatternParametrised::GoldenPatternParametrised(const GoldenPattern* gp): GoldenPattern(gp->key(), gp->getConfig() ) {
  //reset(); called by GoldenPattern()
  meanDistPhi = gp->getMeanDistPhi();
  //std::cout<<key()<<std::endl;
  for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      TLinearFitter linearFitter;
      linearFitter.StoreData(kTRUE);

      std::ostringstream ostr;
      ostr<<"fit_lay_"<<iLayer<<"_refLay_"<<iRefLayer;
      TF1 linearFit("linearFit", "pol3", 0, gp->getPdf()[iLayer][iRefLayer].size());
      linearFitter.SetFormula(linearFit.GetFormula());

      for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
        if(gp->pdfValue(iLayer, iRefLayer, iPdf) > 0) {
          double x = iPdf;
          linearFitter.AddPoint(&x, gp->pdfValue(iLayer, iRefLayer, iPdf));
        }
      }
      int res = linearFitter.Eval();
      //std::cout<<key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" linearFitter result "<<res<<" ";
      if(res == 0) {
        for(unsigned int i = 0; i < pdfParamsAllRef[iLayer][iRefLayer].size(); i++) {
          pdfParamsAllRef[iLayer][iRefLayer][i] = linearFit.GetParameter(i);
          //std::cout<<linearFit.GetParName(i)<<" "<<linearFit.GetParameter(i)<<" ";
        }

      }
      //std::cout<<std::endl;
    }
  }

  generateLuts();

  for(unsigned int iRefLayer=0; iRefLayer < getPdf()[0].size(); ++iRefLayer) {
    for(unsigned int iLayer = 0; iLayer < getPdf().size(); ++iLayer) {
      std::cout<<key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<std::endl;
      for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
        if(gp->pdfValue(iLayer, iRefLayer, iPdf) || pdfValue(iLayer, iRefLayer, iPdf) > 0)
          std::cout<<std::setw(3)<<iPdf<<" org "<<std::setw(3)<<gp->pdfValue(iLayer, iRefLayer, iPdf)<<" param "<<pdfValue(iLayer, iRefLayer, iPdf)<<std::endl;;
      }
    }
  }
}

/*int GoldenPatternParametrised::meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB) const {
  return meanDistPhi[iLayer][iRefLayer][0]; //TODO - choose the right version
  //return ( ( (meanDistPhi[iLayer][iRefLayer][1]*refLayerPhiB)>>myOmtfConfig->nPdfAddrBits() ) + meanDistPhi[iLayer][iRefLayer][0] );
  //assumes that the meanDistPhi[1] is float alpha from the fit to the phiB-phi distribution multiplied by 2^myOmtfConfig->nPdfAddrBits()
}*/

////////////////////////////////////////////////////
////////////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const GoldenPatternParametrised & aPattern){

/*  out <<"GoldenPatternParametrised "<< aPattern.theKey <<std::endl;
  out <<"Number of reference layers: "<<aPattern.meanDistPhi[0].size()
          <<", number of measurement layers: "<<aPattern.pdfAllRef.size()
          <<std::endl;

  if(!aPattern.meanDistPhi.size()) return out;
  if(!aPattern.pdfAllRef.size()) return out;

  out<<"Mean dist phi per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer<<" (";
    for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
      out<<std::setw(3)<<aPattern.meanDistPhi[iLayer][iRefLayer]<<"\t";
    }
    out<<")"<<std::endl;
  }

  if(aPattern.meanDistPhiCounts.size()){
    out<<"Counts number per layer:"<<std::endl;
    for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
      out<<"Ref layer: "<<iRefLayer<<" (";
      for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
        out<<aPattern.meanDistPhiCounts[iLayer][iRefLayer]<<"\t";
      }
      out<<")"<<std::endl;
    }
  }

  unsigned int nPdfAddrBits = 7;
  out<<"PDF per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.pdfAllRef[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer;
    for (unsigned int iLayer=0;iLayer<aPattern.pdfAllRef.size();++iLayer){   
      out<<", measurement layer: "<<iLayer<<std::endl;
      for (unsigned int iRefPhiB=0; iRefPhiB < aPattern.pdfAllRef[iLayer][iRefLayer].size(); ++iRefPhiB) {
        for(unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){
          out<<std::setw(2)<<aPattern.pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf]<<" ";
        }
      }
      out<<std::endl;
    }
  }*/

  return out;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternParametrised::reset(const OMTFConfiguration* omtfConfig){
  GoldenPattern::reset();

  std::vector<floatType> pdf1D(4, 0);
  std::vector<std::vector<floatType> > pdf2D(myOmtfConfig->nRefLayers());
  pdf2D.assign(myOmtfConfig->nRefLayers(), pdf1D);
  pdfParamsAllRef.assign(myOmtfConfig->nLayers(),pdf2D);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
/*int GoldenPatternParametrised::pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB) const {
  if(iBin == 0)
    return 0;

  if((int)iLayer == myOmtfConfig->getRefToLogicNumber()[iRefLayer] && iBin == 64) //FIXME correct 64
      return 63;

  float val = 0;
  unsigned int iBinPow = 1;
  for(unsigned int i = 0; i < pdfParamsAllRef[iLayer][iRefLayer].size(); i++) {
    val += pdfParamsAllRef[iLayer][iRefLayer][i] * iBinPow ;
    iBinPow *= iBin;
  }
  return rint(val);

  //return rint(getA(iLayer, iRefLayer) * (iBin - getB(iLayer, iRefLayer) ) * (iBin - getB(iLayer, iRefLayer) ) + getC(iLayer, iRefLayer) );
}*/

float GoldenPatternParametrised::polynomianValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin) const {
  float val = 0;
  unsigned int iBinPow = 1;
  for(unsigned int i = 0; i < pdfParamsAllRef[iLayer][iRefLayer].size(); i++) {
    val += pdfParamsAllRef[iLayer][iRefLayer][i] * iBinPow ;
    iBinPow *= iBin;
  }
  return val;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternParametrised::generateLuts() {
  for(unsigned int iRefLayer=0; iRefLayer < getPdf()[0].size(); ++iRefLayer) {
    for(unsigned int iLayer = 0; iLayer < getPdf().size(); ++iLayer) {
      setPdfValue(0, iLayer, iRefLayer, 0);
      if((int)iLayer == myOmtfConfig->getRefToLogicNumber()[iRefLayer]) {//FIXME correct 64
          setPdfValue(63, iLayer, iRefLayer, 64);
          continue;
      }
      for(unsigned int iPdf = 1; iPdf < getPdf()[iLayer][iRefLayer].size(); iPdf++) {
        float val = polynomianValue(iLayer, iRefLayer, iPdf);
        if(val < 0)
          val = 0;
        setPdfValue(rint(val), iLayer, iRefLayer, iPdf);
      }

      bool was0 =  false;
      for(unsigned int iPdf = getPdf()[iLayer][iRefLayer].size()/2; iPdf < getPdf()[iLayer][iRefLayer].size(); iPdf++) {
        if(pdfValue(iLayer, iRefLayer, iPdf) <= 0) {
          was0 = true;
        }
        else if(was0) {
          setPdfValue(0, iLayer, iRefLayer, iPdf);
        }
      }

      was0 =  false;
      for(unsigned int iPdf = getPdf()[iLayer][iRefLayer].size()/2; iPdf != 0; iPdf--) {
        if(pdfValue(iLayer, iRefLayer, iPdf) <= 0) {
          was0 = true;
        }
        else if(was0) {
          setPdfValue(0, iLayer, iRefLayer, iPdf);
        }
      }
    }
  }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
int GoldenPatternParametrised::propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer){
  unsigned int iLayer = 2; //MB2
  //if(etaRef>101) iLayer = 7;//RE2
  return phiRef + meanDistPhi[iLayer][iRefLayer][0];
  //FIXME if the meanDistPhiAlpha is non-zero, then meanDistPhi is alone not good for propagation of the phi
  //other value should be used, or the ref_layer phiB should be included
}
