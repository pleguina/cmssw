// RpcKeyToRDumper
//
// Standalone geometry dumper for WP8 follow-up (RPC endcap R-based
// extrapolation, omtf-firmware). Exhaustively iterates every RPC roll in
// the three station/ring=3 endcap chambers OMTF actually reads
// (RE1/3, RE2/3, RE3/3 -- from hwToLogicLayer_0x0210.xml's LayerMap entries
// for logic layers 15-17) and every valid strip in each roll, and for each
// one calls the exact same production code path OMTFinputMaker.cc uses when
// building a real RPC stub (OmtfAngleConverter::getGlobalEtaRpc(rawid,
// strip, r)). Dumps (region,ring,station,sector,subsector,layer,roll,strip,
// etaHw,r) to CSV.
//
// Mirrors CscKeyWGToRDumper.cc's rationale exactly: the omtf-firmware
// phi_extrapolation module only ever has a stub's small discretized eta
// code -- never the real R CMSSW computes from full roll geometry. This
// dump is the exhaustive, geometry-exact source data for deriving an
// analogous (linkID, strip) -> R piecewise table for RPC endcap layers,
// the same way keyWGToR() already does for CSC.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfAngleConverter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"

#include <fstream>
#include <vector>

// checkAndUpdateGeometry() unconditionally calls config->nPhiBins() before
// getGlobalEtaRpc() is ever reached, even though getGlobalEtaRpc() itself
// never reads it -- this dumper only needs the call not to crash, so every
// other pure virtual is stubbed with a value that's never actually
// exercised. Identical stub used by CscKeyWGToRDumper.
class RpcMinimalProcConfig : public ProcConfigurationBase {
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

class RpcKeyToRDumper : public edm::one::EDAnalyzer<> {
public:
  explicit RpcKeyToRDumper(const edm::ParameterSet& cfg)
      : outPath_(cfg.getParameter<std::string>("outPath")),
        cscGeometryEsToken_(esConsumes<CSCGeometry, MuonGeometryRecord>()),
        rpcGeometryEsToken_(esConsumes<RPCGeometry, MuonGeometryRecord>()),
        dtGeometryEsToken_(esConsumes<DTGeometry, MuonGeometryRecord>()) {}

  void analyze(const edm::Event&, const edm::EventSetup& es) override {
    const RPCGeometry& rpcGeometry = es.getData(rpcGeometryEsToken_);

    RpcMinimalProcConfig config;
    MuonGeometryTokens tokens{rpcGeometryEsToken_, cscGeometryEsToken_, dtGeometryEsToken_};
    OmtfAngleConverter angleConverter;
    angleConverter.checkAndUpdateGeometry(es, &config, tokens);

    std::ofstream out(outPath_);
    out << "region,ring,station,sector,subsector,layer,roll,strip,etaHw,r\n";

    for (const RPCRoll* roll : rpcGeometry.rolls()) {
      if (!roll)
        continue;
      RPCDetId id = roll->id();

      // OMTF's hwToLogicLayer_0x0210.xml maps logic layers 15-17 to
      // RE1/3, RE2/3, RE3/3 -- endcap (region != 0), ring 3, station 1-3.
      if (id.region() == 0)
        continue;
      if (id.ring() != 3)
        continue;
      if (id.station() < 1 || id.station() > 3)
        continue;

      const int nstrips = roll->nstrips();
      for (int strip = 1; strip <= nstrips; ++strip) {
        float r = 0;
        int etaHw = angleConverter.getGlobalEtaRpc(id.rawId(), (unsigned int)strip, r);

        out << id.region() << "," << id.ring() << "," << id.station() << "," << id.sector() << "," << id.subsector()
            << "," << id.layer() << "," << id.roll() << "," << strip << "," << etaHw << "," << r << "\n";
      }
    }
    out.close();
    edm::LogInfo("RpcKeyToRDumper") << "Wrote " << outPath_;
  }

private:
  std::string outPath_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryEsToken_;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeometryEsToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryEsToken_;
};

DEFINE_FWK_MODULE(RpcKeyToRDumper);
