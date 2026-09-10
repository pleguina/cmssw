// CscKeyWGToRDumper
//
// Standalone geometry dumper for WP8 (SW-NEW t35 R-based extrapolation
// migration, omtf-firmware). Exhaustively iterates every CSC chamber in the
// four station/ring combinations OMTF actually reads
// (ME1/2, ME1/3, ME2/2, ME3/2 -- from hwToLogicLayer_0x0210.xml's
// LayerMap entries for logic layers 6-9) and every valid wire group in each
// chamber, and for each one calls the exact same production code path
// OMTFinputMaker.cc uses when building a real CSC stub
// (OmtfAngleConverter::getGlobalEta(rawid, digi, r)) with a synthetic
// CSCCorrelatedLCTDigi carrying that wire group. Dumps
// (endcap,station,ring,chamber,keyWG,etaHw,r) to CSV.
//
// This exists to answer a specific firmware question: the omtf-firmware
// phi_extrapolation module only ever has a stub's keyWG (via the same
// coarse etaKeyWG2Code() range-check pattern already used for eta in
// algo/csc_interface/OMTF_csc_interface.cpp) -- never the real R CMSSW
// computes from full chamber geometry. This dump is the exhaustive,
// geometry-exact (not sampled from whichever events happen to appear in
// the verification dataset) source data for deriving an analogous
// keyWG -> R (or keyWG -> R-keyed-extrapolation-factor) piecewise table.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfAngleConverter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"

#include <fstream>
#include <vector>
#include <utility>

// checkAndUpdateGeometry() unconditionally calls config->nPhiBins() before
// getGlobalEta() is ever reached, even though getGlobalEta() itself never
// reads it -- this dumper only needs the call not to crash, so every other
// pure virtual is stubbed with a value that's never actually exercised.
class MinimalProcConfig : public ProcConfigurationBase {
public:
  unsigned int nPhiBins() const override { return 5400; }
  unsigned int nProcessors() const override { return 3; }
  double hwPtToGev(int) const override { return 0; }
  int ptGevToHw(double) const override { return 0; }
  int getProcScalePhi(double, double) const override { return 0; }
  int etaToHwEta(double) const override { return 0; }
  int mb1W2Eta() const override { return 0; }
  int mb2W2Eta() const override { return 0; }
  int mb3W2Eta() const override { return 0; }
  int mb4W2Eta() const override { return 0; }
  unsigned int nLayers() const override { return 18; }
  bool isBendingLayer(unsigned int) const override { return false; }
};

class CscKeyWGToRDumper : public edm::one::EDAnalyzer<> {
public:
  explicit CscKeyWGToRDumper(const edm::ParameterSet& cfg)
      : outPath_(cfg.getParameter<std::string>("outPath")),
        cscGeometryEsToken_(esConsumes<CSCGeometry, MuonGeometryRecord>()),
        rpcGeometryEsToken_(esConsumes<RPCGeometry, MuonGeometryRecord>()),
        dtGeometryEsToken_(esConsumes<DTGeometry, MuonGeometryRecord>()) {}

  void analyze(const edm::Event&, const edm::EventSetup& es) override {
    const CSCGeometry& cscGeometry = es.getData(cscGeometryEsToken_);

    MinimalProcConfig config;
    MuonGeometryTokens tokens{rpcGeometryEsToken_, cscGeometryEsToken_, dtGeometryEsToken_};
    OmtfAngleConverter angleConverter;
    angleConverter.checkAndUpdateGeometry(es, &config, tokens);

    // (endcap, station, ring) pairs OMTF's hwToLogicLayer_0x0210.xml maps to
    // logic layers 6-9 (ME1/3, ME2/2, ME3/2, ME1/2). endcap: 1=+z, 2=-z.
    const std::vector<std::pair<int, int>> stationRings = {{1, 2}, {1, 3}, {2, 2}, {3, 2}};

    std::ofstream out(outPath_);
    out << "endcap,station,ring,chamber,keyWG,etaHw,r\n";

    for (const auto& chamber : cscGeometry.chambers()) {
      CSCDetId chId = chamber->id();
      bool wanted = false;
      for (const auto& sr : stationRings) {
        if (chId.station() == sr.first && chId.ring() == sr.second) {
          wanted = true;
          break;
        }
      }
      if (!wanted)
        continue;

      // Key layer (ALCT/wire-group-bearing layer), matching
      // CSCConstants::KEY_ALCT_LAYER used by getGlobalEta itself.
      const CSCLayer* keyLayer = chamber->layer(3);  // KEY_ALCT_LAYER == 3
      if (!keyLayer)
        continue;
      const CSCLayerGeometry* layerGeom = keyLayer->geometry();
      int nWireGroups = layerGeom->numberOfWireGroups();
      int nStrips = layerGeom->numberOfStrips();
      int midStrip = nStrips / 2;

      for (int wg = 0; wg < nWireGroups; ++wg) {
        CSCCorrelatedLCTDigi digi(/*trknmb*/ 1,
                                   /*valid*/ 1,
                                   /*quality*/ 4,
                                   /*keywire*/ wg,
                                   /*strip*/ midStrip,
                                   /*pattern*/ 10,
                                   /*bend*/ 0,
                                   /*bx*/ 0);
        float r = 0;
        int etaHw = angleConverter.getGlobalEta(chId.rawId(), digi, r);

        out << chId.endcap() << "," << chId.station() << "," << chId.ring() << "," << chId.chamber() << "," << wg
            << "," << etaHw << "," << r << "\n";
      }
    }
    out.close();
    edm::LogInfo("CscKeyWGToRDumper") << "Wrote " << outPath_;
  }

private:
  std::string outPath_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryEsToken_;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeometryEsToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryEsToken_;
};

DEFINE_FWK_MODULE(CscKeyWGToRDumper);
