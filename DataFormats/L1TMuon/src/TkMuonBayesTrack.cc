/*
 * TkMuonBayesTrack.cc
 *
 *  Created on: Mar 15, 2019
 *      Author: Author: Karol Bunkowski kbunkow@cern.ch
 */

#include "DataFormats/L1TMuon/interface/TkMuonBayesTrack.h"

namespace l1t {
  TkMuonBayesTrack::TkMuonBayesTrack() : L1Candidate() {}

  TkMuonBayesTrack::TkMuonBayesTrack(const LorentzVector& p4) : L1Candidate(p4) {}

  TkMuonBayesTrack::TkMuonBayesTrack(const edm::Ptr<L1TTTrackType> ttTrackPtr)
      : L1Candidate(LorentzVector(ttTrackPtr->momentum().x(),
                                  ttTrackPtr->momentum().y(),
                                  ttTrackPtr->momentum().z(),
                                  ttTrackPtr->momentum().mag())),
        ttTrackPtr(ttTrackPtr) {}

  TkMuonBayesTrack::~TkMuonBayesTrack() {}

}  //end of namespace l1t
