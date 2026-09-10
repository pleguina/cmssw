// RpcKeyToRDumper
//
// Standalone geometry dumper for WP8 follow-up (RPC endcap R-based
// extrapolation, omtf-firmware). Exhaustively iterates every RPC roll in
// the three station/ring=3 endcap chambers OMTF actually reads
// (RE1/3, RE2/3, RE3/3 -- from hwToLogicLayer_0x0210.xml's LayerMap entries
// for logic layers 15-17) and every valid strip in each roll, and for each
// one calls the exact same production code path OMTFinputMaker.cc uses when
// building a real RPC stub: OmtfAngleConverter::getGlobalEtaRpc(rawid,
// strip, r) for r/eta, and OmtfAngleConverter::getProcessorPhi(phiZero=0,
// part, rollId, strip, strip) for phi (dumped with phiZero=0 so it carries
// the raw, processor-independent angle bin -- OMTFinputMaker.cc's own
// getProcessorPhiZero(iProcessor) = nPhiBins/nProcessors*iProcessor +
// nPhiBins/24 is a closed-form per-processor constant, cheap to apply in
// firmware directly rather than needing per-processor dump columns).
// Dumps (region,ring,station,sector,subsector,layer,roll,strip,etaHw,r,
// phiBin) to CSV.
//
// Mirrors CscKeyWGToRDumper.cc's rationale exactly: the omtf-firmware
// phi_extrapolation module only ever has a stub's small discretized eta
// code -- never the real R CMSSW computes from full roll geometry. This
// dump is the exhaustive, geometry-exact source data for deriving an
// analogous (linkID, strip) -> R/eta/phi piecewise table for RPC endcap
// layers, the same way keyWGToR() already does for CSC.

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
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfAngleConverter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"

#include <cmath>
#include <fstream>
#include <vector>

// checkAndUpdateGeometry() unconditionally calls config->nPhiBins() before
// getGlobalEtaRpc() is ever reached, even though getGlobalEtaRpc() itself
// never reads it -- this dumper only needs the call not to crash, so every
// other pure virtual is stubbed with a value that's never actually
// exercised. Identical stub used by CscKeyWGToRDumper, EXCEPT
// etaToHwEta(): implemented for real here (abs(eta/etaUnit_), the exact
// formula OmtfAngleConverter::getGlobalEtaRpc uses for
// StubEtaEncoding::valueP1Scale/valueP2Scale, per OmtfAngleConverter.cc)
// so this dump captures real SW-NEW-correct eta, not the legacy 'bits'
// one-hot mask CscKeyWGToRDumper's identical-looking stub falls back to
// when etaToHwEta returns 0 (that mistake was caught and NOT reused here
// -- see reference_manifest.yaml's WP8 entries).
class RpcMinimalProcConfig : public ProcConfigurationBase {
public:
  unsigned int nPhiBins() const override { return 5400; }
  unsigned int nProcessors() const override { return 3; }
  double hwPtToGev(int) const override { return 0; }
  int ptGevToHw(double) const override { return 0; }
  int getProcScalePhi(double, double) const override { return 0; }
  int etaToHwEta(double eta) const override { return (int)std::lround(std::abs(eta) / etaUnit_); }
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
    config.setStubEtaEncoding(ProcConfigurationBase::StubEtaEncoding::valueP1Scale);
    MuonGeometryTokens tokens{rpcGeometryEsToken_, cscGeometryEsToken_, dtGeometryEsToken_};
    OmtfAngleConverter angleConverter;
    angleConverter.checkAndUpdateGeometry(es, &config, tokens);

    std::ofstream out(outPath_);
    out << "region,ring,station,sector,subsector,layer,roll,strip,etaHw,r,phiBin\n";

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
        // phiZero=0: raw, processor-independent angle bin (see file header
        // -- the real per-processor phiZero is a cheap closed-form
        // subtraction applied in firmware, not baked into this dump).
        int phiBin = angleConverter.getProcessorPhi(
            /*phiZero=*/0, l1t::tftype::omtf_pos, id, (unsigned int)strip, (unsigned int)strip);

        out << id.region() << "," << id.ring() << "," << id.station() << "," << id.sector() << "," << id.subsector()
            << "," << id.layer() << "," << id.roll() << "," << strip << "," << etaHw << "," << r << "," << phiBin
            << "\n";
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
