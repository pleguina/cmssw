/*
 * DataDumper.h
 *
 *  Created on: Nov 5, 2021
 *      Author: kbunkow
 */

#ifndef INTERFACE_DATADUMPER_H_
#define INTERFACE_DATADUMPER_H_

#include "PreTrackMatchedMuon.h"
#include "ConvertedTTTrack.h"
#include "MuonROI.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"

#include "TTree.h"
#include "TFile.h"

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/archive/xml_oarchive.hpp>


namespace l1t {

template<class Archive>
void serialize(Archive & ar, MuonStub& stub, const unsigned int version) {
  //only saves the data, will not load them - but this is not needed
  auto etaRegion = stub.etaRegion();
  ar & BOOST_SERIALIZATION_NVP(etaRegion);

  auto phiRegion = stub.phiRegion();
  ar & BOOST_SERIALIZATION_NVP(phiRegion);

  auto depthRegion = stub.depthRegion();
  ar & BOOST_SERIALIZATION_NVP(depthRegion);

  auto coord1 = stub.coord1();
  ar & BOOST_SERIALIZATION_NVP(coord1);

  auto coord2 = stub.coord2();
  ar & BOOST_SERIALIZATION_NVP(coord2);

  auto id = stub.id();
  ar & BOOST_SERIALIZATION_NVP(id);

  auto quality = stub.quality();
  ar & BOOST_SERIALIZATION_NVP(quality);

  auto bxNum = stub.bxNum();
  ar & BOOST_SERIALIZATION_NVP(bxNum);

  auto eta1 = stub.eta1();
  ar & BOOST_SERIALIZATION_NVP(eta1);

  auto eta2 = stub.eta2();
  ar & BOOST_SERIALIZATION_NVP(eta2);

  auto etaQuality = stub.etaQuality();
  ar & BOOST_SERIALIZATION_NVP(etaQuality);

  auto type = stub.type();
  ar & BOOST_SERIALIZATION_NVP(type);
}

};

BOOST_CLASS_VERSION(l1t::MuonStub, 1);

namespace Phase2L1GMT {

template<class Archive>
void serialize(Archive & ar, ConvertedTTTrack& mu, const unsigned int version) {
  auto charge = mu.charge();
  ar & BOOST_SERIALIZATION_NVP(charge);

  auto curvature = mu.curvature();
  ar & BOOST_SERIALIZATION_NVP(curvature);

  auto abseta = mu.abseta();
  ar & BOOST_SERIALIZATION_NVP(abseta);

  auto pt = mu.pt();
  ar & BOOST_SERIALIZATION_NVP(pt);

  auto eta = mu.eta();
  ar & BOOST_SERIALIZATION_NVP(eta);

  auto phi = mu.phi();
  ar & BOOST_SERIALIZATION_NVP(phi);

  auto z0 = mu.z0();
  ar & BOOST_SERIALIZATION_NVP(z0);

  auto d0 = mu.d0();
  ar & BOOST_SERIALIZATION_NVP(d0);

  auto quality = mu.quality();
  ar & BOOST_SERIALIZATION_NVP(quality);

/*  auto offline_pt = mu.offline_pt();
  ar & BOOST_SERIALIZATION_NVP(offline_pt);

  auto offline_eta = mu.offline_eta();
  ar & BOOST_SERIALIZATION_NVP(offline_eta);

  auto offline_phi = mu.offline_phi();
  ar & BOOST_SERIALIZATION_NVP(offline_phi);*/
}

template<class Archive>
void serialize(Archive & ar, PreTrackMatchedMuon& mu, const unsigned int version) {

  //only saves the data, will not load them - but this is not needed
  auto curvature = mu.curvature();
  ar & BOOST_SERIALIZATION_NVP(curvature);

  auto charge = mu.charge();
  ar & BOOST_SERIALIZATION_NVP(charge);

  auto pt = mu.pt();
  ar & BOOST_SERIALIZATION_NVP(pt);

  auto eta = mu.eta();
  ar & BOOST_SERIALIZATION_NVP(eta);

  auto phi = mu.phi();
  ar & BOOST_SERIALIZATION_NVP(phi);

  auto z0 = mu.z0();
  ar & BOOST_SERIALIZATION_NVP(z0);

  auto d0 = mu.d0();
  ar & BOOST_SERIALIZATION_NVP(d0);

  auto beta = mu.beta();
  ar & BOOST_SERIALIZATION_NVP(beta);

  auto isGlobalMuon = mu.isGlobalMuon();
  ar & BOOST_SERIALIZATION_NVP(isGlobalMuon);

  auto quality = mu.quality();
  ar & BOOST_SERIALIZATION_NVP(quality);

/*  auto offline_pt = mu.offline_pt();
  ar & BOOST_SERIALIZATION_NVP(offline_pt);

  auto offline_eta = mu.offline_eta();
  ar & BOOST_SERIALIZATION_NVP(offline_eta);

  auto offline_phi = mu.offline_phi();
  ar & BOOST_SERIALIZATION_NVP(offline_phi);*/

  auto stubID0 = mu.stubID0();
  ar & BOOST_SERIALIZATION_NVP(stubID0);

  auto stubID1 = mu.stubID1();
  ar & BOOST_SERIALIZATION_NVP(stubID1);

  auto stubID2 = mu.stubID2();
  ar & BOOST_SERIALIZATION_NVP(stubID2);

  auto stubID3 = mu.stubID3();
  ar & BOOST_SERIALIZATION_NVP(stubID3);

  auto stubID4 = mu.stubID4();
  ar & BOOST_SERIALIZATION_NVP(stubID4);
}

struct TrackMatchedMuonRecord {
  const unsigned int TF_LAYER_CNT = 5;

  //from trackingParticle
  unsigned int eventNum = 0;
  short tpEvent = -1;
  float tpPt = 0, tpEta = 0, tpPhi = 0, tpBeta = 0;
  short tpCharge = 0;
  int tpType = 0;

  //0 - no match, 1 - very loose, 2 - loose , 3
  char matching = 0;

  //from ttTrack
  short tttCurvature = 0;
  char tttCharge = 0;
  unsigned short tttPt = 0;
  short tttEta = 0;
  short tttPhi = 0;
  short tttZ0 = 0;
  short tttD0 = 0;

  //from TTTrack_TrackWord
  char tttChi2rphi = 0; //4 bits
  char tttChi2rz = 0; //4 bits
  char tttBendChi2 = 0; //3 bits
  char tttQualityMVA = 0; //3 bits
  char tttOtherMVA = 0; //6 bits

  //from GMT
  uint gmtBeta = 0;
  bool isGlobal = 0;
  uint quality = 0;

  unsigned short hitsValid = 0;

  std::vector<unsigned char> deltaCoords1;
  std::vector<unsigned char> deltaCoords2;

  std::vector<unsigned char> deltaEtas1;
  std::vector<unsigned char> deltaEtas2;

  std::vector<char> stubTiming;
  std::vector<char> stubType;


  //std::vector<unsigned char> stubTiming2;

  TrackMatchedMuonRecord():
    deltaCoords1(TF_LAYER_CNT), deltaCoords2(TF_LAYER_CNT),
    deltaEtas1(TF_LAYER_CNT), deltaEtas2(TF_LAYER_CNT),
    stubTiming(TF_LAYER_CNT), stubType(TF_LAYER_CNT) {}
};


class PreTrackMatchedMuonProcessor {
public:
  PreTrackMatchedMuonProcessor() {};
  virtual ~PreTrackMatchedMuonProcessor() {};

  virtual void eventBegin() = 0;

  virtual void collectNonant(int nonant, std::vector<MuonROI>& rois,
      std::vector<ConvertedTTTrack>& convertedTracks,
      std::vector<PreTrackMatchedMuon>& muons,
      std::vector<PreTrackMatchedMuon>& muCleaned) = 0;

  virtual void process(PreTrackMatchedMuon& preTrackMatchedMuon) = 0;

  virtual void writeToXml() = 0;
};

/*
 * if dumpToRoot == true the data are dumped to root file opened by the the edm::Service<TFileService> fs,
 * so it must be opened in the python
 * if dumpToXml == true it produces xml files that are useful for testing the firmware
 */
class DataDumper: public PreTrackMatchedMuonProcessor {
public:
  DataDumper(const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& ttTrackMCTruthToken,
      const edm::EDGetTokenT< std::vector< TrackingParticle > >& trackingParticleToken, bool dumpToRoot, bool dumpToXml);

  virtual ~DataDumper();

  void getHandles(const edm::Event& event);

  void eventBegin() override {
    evntCnt++;
    eventXml.clear();
  }

  void process(PreTrackMatchedMuon& preTrackMatchedMuon) override;

  void collectNonant(int nonant, std::vector<MuonROI>& rois,
      std::vector<ConvertedTTTrack>& convertedTracks,
      std::vector<PreTrackMatchedMuon>& muons,
      std::vector<PreTrackMatchedMuon>& muCleaned) override;

  void writeToXml() override;

private:
  void initializeTTree();

  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken;
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > mcTruthTTTrackHandle;

  edm::EDGetTokenT< std::vector< TrackingParticle > > trackingParticleToken;
  edm::Handle< std::vector< TrackingParticle > > trackingParticleHandle;

  bool dumpToRoot = false;
  bool dumpToXml = false;

  std::vector<edm::Ptr< TrackingParticle > > muonTrackingParticles;
  bool muonTrackingParticlesFilled = false;

  TTree* rootTree = nullptr;

  TrackMatchedMuonRecord record;

  unsigned int evntCnt = 0;

  //TH1I* ptGenPos = nullptr;
  //TH1I* ptGenNeg = nullptr;

  class NonantXml {
  public:
    NonantXml() {
    }

    std::vector<std::vector<l1t::MuonStub> > rois; //[[roiIndex][stubIndex]

    std::vector<ConvertedTTTrack> convertedTracks;

    std::vector<PreTrackMatchedMuon> muons;

    std::vector<PreTrackMatchedMuon> muCleaned;

  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(rois);
      ar & BOOST_SERIALIZATION_NVP(convertedTracks);
      ar & BOOST_SERIALIZATION_NVP(muons);
      ar & BOOST_SERIALIZATION_NVP(muCleaned);
    }
  };

  class EventXml {
  public:
    EventXml() : nonants(9) {

    }

    std::vector<NonantXml > nonants; //[nonant]

    void clear() {
      nonants.clear();
      nonants.resize(9);
      empty = true;
    }

    bool empty = true;
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(nonants);
    }
  };

  EventXml eventXml;
};


} /* namespace Phase2L1GMT */


#endif /* INTERFACE_DATADUMPER_H_ */
