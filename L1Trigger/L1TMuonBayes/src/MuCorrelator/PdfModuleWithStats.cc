/*
 * PdfModuleWithStats.cc
 *
 *  Created on: Feb 4, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/PdfModuleWithStats.h"

PdfModuleWithStats::PdfModuleWithStats(MuCorrelatorConfigPtr& config): PdfModule(config), pdfHists(config->nLayers() ) {
  for(unsigned int iLayer = 0; iLayer < coefficients.size(); ++iLayer) {
    for(unsigned int iEtaBin = 0; iEtaBin < coefficients[iLayer].size(); ++iEtaBin) {
      pdfHists[iLayer].emplace_back();
      for(unsigned int iRefLayer = 0; iRefLayer < coefficients[iLayer][iEtaBin].size(); ++iRefLayer) {
        std::ostringstream name;
                name<<"pdfHist_layer_"<<iLayer<<"_eta_"<<iEtaBin<<"_refLayer_"<<iRefLayer;
        pdfHists[iLayer][iEtaBin].emplace_back(std::make_unique<TH2I>(name.str(). c_str(), name.str(). c_str(),
                                                           config->nPtBins(), 0, config->nPtBins(), 1000, 0, 1000));
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

void PdfModuleWithStats::write() const {
  for(unsigned int iLayer = 0; iLayer < pdfHists.size(); ++iLayer) {
    for(unsigned int iEtaBin = 0; iEtaBin < pdfHists[iLayer].size(); ++iEtaBin) {
      for(unsigned int iRefLayer = 0; iRefLayer < pdfHists[iLayer][iEtaBin].size(); ++iRefLayer) {
        pdfHists.at(iLayer).at(iEtaBin).at(iRefLayer)->Write();
      }
    }
  }
}
