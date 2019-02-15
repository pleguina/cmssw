/*
 * PdfModuleWithStats.cc
 *
 *  Created on: Feb 4, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/PdfModuleWithStats.h"

#include "TFitResultPtr.h"
#include "TFitResult.h"

PdfModuleWithStats::PdfModuleWithStats(MuCorrelatorConfigPtr& config): PdfModule(config), pdfHists(config->nLayers() ) {
  TFileDirectory subDir = fs->mkdir("pdfs");

  for(unsigned int iLayer = 0; iLayer < coefficients.size(); ++iLayer) {
    for(unsigned int iEtaBin = 0; iEtaBin < coefficients[iLayer].size(); ++iEtaBin) {
      pdfHists[iLayer].emplace_back();
      for(unsigned int iRefLayer = 0; iRefLayer < coefficients[iLayer][iEtaBin].size(); ++iRefLayer) {
        std::ostringstream name;
                name<<"pdfHist_layer_"<<iLayer<<"_eta_"<<iEtaBin<<"_refLayer_"<<iRefLayer;
        pdfHists[iLayer][iEtaBin].emplace_back(subDir.make<TH2I>(name.str(). c_str(), name.str(). c_str(),
                                                           config->nPtBins(), 0, config->nPtBins(), 1100, -100, 1000));
      }
    }

    //[layer][etaBin][refLayer](ptBin, pdfBin)
  }
}

PdfModuleWithStats::~PdfModuleWithStats() {
  // TODO Auto-generated destructor stub
}

float PdfModuleWithStats::getPdfVal(unsigned int layer, unsigned int etaBin, unsigned int refLayer, unsigned int ptBin, int pdfBin) {
  pdfHists.at(layer).at(etaBin).at(refLayer)->Fill(ptBin, pdfBin);
  return PdfModule::getPdfVal(layer, etaBin, refLayer, ptBin, pdfBin);
}

/* whith file service it is not needed
void PdfModuleWithStats::write() const {
  for(unsigned int iLayer = 0; iLayer < pdfHists.size(); ++iLayer) {
    for(unsigned int iEtaBin = 0; iEtaBin < pdfHists[iLayer].size(); ++iEtaBin) {
      for(unsigned int iRefLayer = 0; iRefLayer < pdfHists[iLayer][iEtaBin].size(); ++iRefLayer) {
        pdfHists.at(iLayer).at(iEtaBin).at(iRefLayer)->Write();
      }
    }
  }
}*/


void PdfModuleWithStats::generateCoefficients() const {
  //fs->cd();
  TFileDirectory subDir = fs->mkdir("pdfsProj"); //does not work, why?
  //subDir.cd();

  for(unsigned int iLayer = 0; iLayer < pdfHists.size(); ++iLayer) {
    for(unsigned int iEtaBin = 0; iEtaBin < pdfHists[iLayer].size(); ++iEtaBin) {
      for(unsigned int iRefLayer = 0; iRefLayer < pdfHists[iLayer][iEtaBin].size(); ++iRefLayer) {
        auto pdfHist = pdfHists.at(iLayer).at(iEtaBin).at(iRefLayer);
        for(int ptBin = 1; ptBin < pdfHist->GetXaxis()->GetNbins(); ++ptBin) {
          //Normalize pdf in each ptBin separately, to get p(pdBin | pt, eta)
          std::ostringstream ostr;
          ostr<<pdfHist->GetName()<<"_ptBin_"<<ptBin;
          TH1D* pdfHistInPtBin = pdfHist->ProjectionY(ostr.str().c_str(), ptBin, ptBin);

          pdfHistInPtBin->Sumw2();
          if(pdfHistInPtBin->Integral() <= 0) {
            std::cout<<__FUNCTION__<<": "<<__LINE__<<" iLayer "<<iLayer<<" iEtaBin "<<iEtaBin<<" iRefLayer "<<iRefLayer<<" ptBin "<<ptBin<<" - no entries, coefficients not calculated"<<std::endl;
            continue;
          }

          pdfHistInPtBin->Scale(1./pdfHistInPtBin->Integral());

          //pdfHistInPtBin->Write(); //TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,,
          //          continue;

          for(int iBinPdf = 0; iBinPdf < pdfHistInPtBin->GetXaxis()->GetNbins(); iBinPdf++) {
            double pdfVal = pdfHistInPtBin->GetBinContent(iBinPdf);

            const double minPlog =  log(config->minPdfVal());
            const double pdfMaxLogVal = config->pdfMaxLogValue();

            double logPdf = 0;
            double error = 0;
            if(pdfVal >= config->minPdfVal()/2.)  //10 is for removing points with small statistics - TODO tune
            //if(pdfVal > 0)
            {
              logPdf = pdfMaxLogVal - log(pdfVal) / minPlog * pdfMaxLogVal;

              error = 1./pdfVal / minPlog * pdfMaxLogVal * pdfHistInPtBin->GetBinError(iBinPdf);
            }

            pdfHistInPtBin->SetBinContent(iBinPdf, logPdf);
            pdfHistInPtBin->SetBinError(iBinPdf, error);
          }

          //TGraph* gr = new TGraph(iPoint, x, y);

          try {
            TFitResultPtr fitResult = pdfHistInPtBin->Fit("pol2","S");
            Int_t fitStatus = fitResult;
            if(fitStatus == 0) {
              //fitResult->Clear(); is this needed?

              double a0 = fitResult->Value(0);
              double a1 = fitResult->Value(1);
              double a2 = fitResult->Value(2);

              std::cout<<__FUNCTION__<<": "<<__LINE__<<" iLayer "<<iLayer<<" iEtaBin "<<iEtaBin<<" iRefLayer "<<iRefLayer<<" ptBin "<<ptBin
                  <<" a0 "<<a0<<" a1 "<<a1<<" a2 "<<a2<<std::endl;

            }
            else {
              std::cout<<__FUNCTION__<<": "<<__LINE__<<" iLayer "<<iLayer<<" iEtaBin "<<iEtaBin<<" iRefLayer "<<iRefLayer<<" ptBin "<<ptBin
                                <<" - fit failed!!!!!!!!!!!!!!!!!! "<<std::endl;
            }
          }
          catch (std::exception& err) {
            std::cout<<__FUNCTION__<<": "<<__LINE__<<" iLayer "<<iLayer<<" iEtaBin "<<iEtaBin<<" iRefLayer "<<iRefLayer<<" ptBin "<<ptBin
                                 <<" - fit failed "<<err.what()<<std::endl;
          }

          pdfHistInPtBin->Write();
        }
      }
    }
  }
}
