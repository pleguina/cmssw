/*
 * PdfModuleWithStats.h
 *
 *  Created on: Feb 4, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef INTERFACE_MUCORRELATOR_PDFMODULEWITHSTATS_H_
#define INTERFACE_MUCORRELATOR_PDFMODULEWITHSTATS_H_

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/PdfModule.h"

#include "TH2I.h"

class PdfModuleWithStats: public PdfModule {
public:
  PdfModuleWithStats(MuCorrelatorConfigPtr& config);

  virtual ~PdfModuleWithStats();

  virtual float getPdfVal(unsigned int layer, unsigned int etaBin, unsigned int refLayer, unsigned int ptBin, int pdfBin);

  //writes pdfHists to current root file
  virtual void write() const;
private:
  //[layer][etaBin][refLayer](ptBin, pdfBin)
  std::vector<std::vector<std::vector< std::unique_ptr<TH2I> > > > pdfHists;
};

#endif /* INTERFACE_MUCORRELATOR_PDFMODULEWITHSTATS_H_ */
