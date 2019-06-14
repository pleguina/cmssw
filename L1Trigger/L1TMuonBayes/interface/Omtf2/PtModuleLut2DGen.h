/*
 * PtModuleLut2DGen.h
 *
 *  Created on: May 22, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF2_PTMODULELUT2DGEN_H_
#define INTERFACE_OMTF2_PTMODULELUT2DGEN_H_

#include "L1Trigger/L1TMuonBayes/interface/Omtf2/PtModuleLut2D.h"

//#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class PtModuleLut2DGen: public PtModuleLut2D {
public:
  PtModuleLut2DGen(const ProcConfigurationBase* config);

  virtual ~PtModuleLut2DGen();

  enum Mode {
/*    calcualtePtMeans,
    calcualtePtStdDevs,
    calcualteDispMeans,
    calcualteDispStdDevs*/
    calcualteMeans,
    calcualteStdDevs
  };

  enum SampleType {
    muons,
    displacedMuons
  };

  void prepare(Mode mode, SampleType sampleType);

  //the simulated muon corresponding to the muon stubs should be given for each event, assuming that there is only one in the event
/*  void setSimMuon(const SimTrack* simMuon) {
    this->simMuon = simMuon;
  }*/

  //the simulated muon corresponding to the muon stubs should be given for each event, assuming that there is only one in the event
  void setGenMuon(const reco::GenParticle* genMuon) {
    this->genMuon = genMuon;
  }

  //empty muon is returned
  virtual AlgoMuon2Ptr processStubs(const MuonStubsInput& muonStubs);

  virtual void endJob();

  void fillEmptyCells(unsigned int iOut, unsigned int iOutLike);

  void calcualteLogLikelihoods(unsigned int iOut);

private:
  //const SimTrack* simMuon = nullptr;
  const reco::GenParticle* genMuon = nullptr;

  Mode mode = Mode::calcualteMeans;
  SampleType sampleType = SampleType::muons;

};

#endif /* INTERFACE_OMTF2_PTMODULELUT2DGEN_H_ */
