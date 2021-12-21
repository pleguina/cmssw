#ifndef PHASE2GMT_CONEVRTEDTTRACK
#define PHASE2GMT_CONEVRTEDTTRACK
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace Phase2L1GMT {

  class ConvertedTTTrack {
  public:
    ConvertedTTTrack(const uint& charge,
                     const int& curvature,
                     const uint& abseta,
                     const uint& pt,
                     const int& eta,
                     const int& phi,
                     const int& z0,
                     const int& d0,
                     const int& quality,
                     const ap_uint<96>& word)
        : charge_(charge),
          curvature_(curvature),
          abseta_(abseta),
          pt_(pt),
          eta_(eta),
          phi_(phi),
          z0_(z0),
          d0_(d0),
          quality_(quality),
          word_(word) {}

    const uint charge() const { return charge_; }

    const int curvature() const { return curvature_; }
    const uint abseta() const { return abseta_; }

    const uint pt() const { return pt_; }

    const int eta() const { return eta_; }
    const int phi() const { return phi_; }

    void setPhi(int phi) { phi_ = phi; }

    const int z0() const { return z0_; }
    const int d0() const { return d0_; }
    const int quality() const { return quality_; }
    const float offline_pt() const { return offline_pt_; }
    const float offline_eta() const { return offline_eta_; }
    const float offline_phi() const { return offline_phi_; }

    const ap_uint<96>& word() const { return word_; }
    void setOfflineQuantities(float pt, float eta, float phi) {
      offline_pt_ = pt;
      offline_eta_ = eta;
      offline_phi_ = phi;
    }

    friend std::ostream& operator<<(std::ostream& out, const ConvertedTTTrack& track) {
      out << "converted track : charge " << track.charge_ << " curvature " << track.curvature_
          << " pt=" << track.offline_pt_ << " GeV = " << track.pt_ << " eta " << track.offline_eta_ << " = " << track.eta_
          << " phi=" << track.offline_phi_ << " = " << track.phi_ << " z0 " << track.z0_ << " d0 " << track.d0_
          << " quality " << track.quality_<<" abseta_ "<<track.abseta_;

      return out;
    }

    void printWord() const {
      LogTrace("ConvertedTTTrack") << "converted track : word=" << std::setfill('0') << std::setw(8) << std::hex
                                   << (long long unsigned int)((word_ >> 64).to_uint64()) << std::setfill('0')
                                   << std::setw(16) << std::hex
                                   << (long long unsigned int)((word_ & 0xffffffffffffffff).to_uint64());
    }

    void setTrkPtr(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& trkPtr) { trkPtr_ = trkPtr; }

    const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr() const { return trkPtr_; }

  private:
    uint charge_;
    int curvature_;
    uint abseta_;
    uint pt_;
    int eta_;
    int phi_;
    int z0_;
    int d0_;
    uint quality_;
    float offline_pt_ = 0;
    float offline_eta_ = 0;
    float offline_phi_ = 0;
    ap_uint<96> word_;

    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr_;
  };
}  // namespace Phase2L1GMT

#endif
