/*
 * TTAlgoMuon.cc
 *
 *  Created on: Oct 19, 2018
 *      Author: kbunkow
 */


#include "L1Trigger/L1TMuonOverlap/interface/TTAlgoMuon.h"

TrackingTriggerTrack::TrackingTriggerTrack(const SimTrack& simMuon) {
  eta = simMuon.momentum().eta();
  phi = simMuon.momentum().phi();
  pt = simMuon.momentum().pt();
  charge = simMuon.type()/-13;
}

