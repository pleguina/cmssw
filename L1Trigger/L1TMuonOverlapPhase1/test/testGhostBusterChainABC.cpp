//============================================================================
// Standalone driver for GhostBuster::select (the SW reference GhostBuster
// implementation, GhostBuster.cc, instantiated by default from
// OMTFProcessor::ghostBust) exercising the "chain A-B-C" directed test case
// required by omtf-firmware's docs/OMTF_VERIFICATION_OPTIMIZED_PLAN.md
// Sec.8.1: A close to B, B close to C, A NOT close to C.
//
// GhostBuster.cc itself carries an unresolved TODO left by its author
// (initials "KB") asking exactly this question: the SW's suppression is
// greedy-sequential (a candidate is only ever compared against the
// already-*accepted* list, so a killed candidate can never itself kill a
// third one) -- does the firmware (algo/ghostBuster/OMTF_ghostBuster.cpp,
// killGhosts(), a full O(N^2) pairwise matrix over the *original*
// candidates) behave the same way? This binary answers that question for
// the SW side; the matching HW-side scenario 2 lives in
// omtf-firmware's arc/verify/tests/tb_ghostbuster.cpp. Both use the same
// phi/rank values so the two runs are directly comparable -- see
// docs/OMTF_VERIFICATION_GHOSTBUSTER_CHAIN_ABC.md in omtf-firmware for the
// side-by-side result.
//
// Build & run (from a CMSSW_14_2_0_pre2 dev area with this package checked
// out under src/L1Trigger/L1TMuonOverlapPhase1):
//   scram b
//   testGhostBusterChainABC
//============================================================================
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GhostBuster.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/AlgoMuon.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPattern.h"

#include <cassert>
#include <iostream>

namespace {

// GhostBuster::select copies each surviving AlgoMuon (`new AlgoMuon(**it1)`)
// but only shallow-copies the goldenPaternConstr pointer, so the GoldenPattern
// backing each candidate must outlive the whole test -- intentionally never
// deleted here, this is a short-lived one-shot binary, not a leak to fix.
AlgoMuonPtr makeCandidate(int procPhi, unsigned int firedLayerCnt, unsigned int patternNumber,
                           unsigned int refHitNumber) {
  GoldenPatternResult gpr;
  gpr.setFiredLayerCnt(firedLayerCnt);
  gpr.setPhi(procPhi);
  // pt=10 (>0) is all isValid()/getPtConstr() need; the actual pt value is
  // irrelevant to GhostBuster::select, which ranks by Q/Disc/PatternNumConstr/
  // RefHitNumber and vetoes by phi only.
  auto* gp = new GoldenPattern(Key(/*eta*/ 0, /*pt*/ 10, /*charge*/ 1, patternNumber),
                                /*nLayers*/ 1, /*nRefLayers*/ 1, /*nPdfAddrBits*/ 1);
  return std::make_shared<AlgoMuon>(gpr, gp, refHitNumber);
}

void printCand(const char* label, const AlgoMuonPtr& c) {
  std::cout << label << ": Q=" << c->getQ() << " phi=" << c->getPhi() << " valid=" << c->isValid() << "\n";
}

}  // namespace

int main() {
  OMTFConfiguration omtfConfig;  // procPhiToGmtPhi() is a pure function of its argument; configure() not needed.

  // procPhi values chosen so that, after procPhiToGmtPhi() conversion (the
  // same GMT phi scale killGhosts() operates on in hardware phi units),
  // gmtPhi(A)=100, gmtPhi(B)=105, gmtPhi(C)=111 -- |A-B|=5<8, |B-C|=6<8,
  // |A-C|=11>=8: A-B and B-C are ghost pairs, A-C is not. Same numbers used
  // in tb_ghostbuster.cpp's scenario 2 for a direct side-by-side comparison.
  AlgoMuonPtr A = makeCandidate(/*procPhi*/ 940, /*firedLayerCnt(Q)*/ 5, /*patternNumber*/ 1, /*refHitNumber*/ 0);
  AlgoMuonPtr B = makeCandidate(/*procPhi*/ 988, /*firedLayerCnt(Q)*/ 4, /*patternNumber*/ 2, /*refHitNumber*/ 1);
  AlgoMuonPtr C = makeCandidate(/*procPhi*/ 1045, /*firedLayerCnt(Q)*/ 3, /*patternNumber*/ 3, /*refHitNumber*/ 2);

  std::cout << "gmtPhi: A=" << omtfConfig.procPhiToGmtPhi(A->getPhi())
            << " B=" << omtfConfig.procPhiToGmtPhi(B->getPhi())
            << " C=" << omtfConfig.procPhiToGmtPhi(C->getPhi()) << "\n";
  std::cout << "|A-B|=" << std::abs(omtfConfig.procPhiToGmtPhi(A->getPhi()) - omtfConfig.procPhiToGmtPhi(B->getPhi()))
            << " |B-C|=" << std::abs(omtfConfig.procPhiToGmtPhi(B->getPhi()) - omtfConfig.procPhiToGmtPhi(C->getPhi()))
            << " |A-C|=" << std::abs(omtfConfig.procPhiToGmtPhi(A->getPhi()) - omtfConfig.procPhiToGmtPhi(C->getPhi()))
            << " (veto window: <8)\n";

  // Deliberately out of rank order in the input vector -- select() sorts by
  // goodness itself (GhostBuster.cc:35), so this also exercises that sort.
  AlgoMuons refHitCands{C, A, B};

  GhostBuster gb(&omtfConfig);
  AlgoMuons out = gb.select(refHitCands);

  std::cout << "-- GhostBuster::select output (" << out.size() << " candidates) --\n";
  for (size_t i = 0; i < out.size(); i++) {
    std::cout << "out[" << i << "]: ";
    printCand("", out[i]);
  }

  assert(out.size() == 3 && "select() always pads/caps refHitCleanCands to exactly 3");
  std::cout << "\nRESULT: SW ghost-buster keeps "
            << ((out[0]->getQ() > 0) + (out[1]->getQ() > 0) + (out[2]->getQ() > 0)) << " live candidate(s) of A,B,C.\n";

  return 0;
}
