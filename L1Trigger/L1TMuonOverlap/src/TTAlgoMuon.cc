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

TrackingTriggerTrack::TrackingTriggerTrack(const TTTrack< Ref_Phase2TrackerDigi_>& ttTRack, int l1Tk_nPar) {
    pt   = ttTRack.getMomentum(l1Tk_nPar).perp();
    eta  = ttTRack.getMomentum(l1Tk_nPar).eta();
    phi  = ttTRack.getMomentum(l1Tk_nPar).phi();

    charge = (ttTRack.getRInv() > 0 ? 1 : -1); //ttTRack.ge //where is the charge???? TODO
}
