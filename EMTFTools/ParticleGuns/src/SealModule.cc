#include "FWCore/Framework/interface/MakerMacros.h"

#include "EMTFTools/ParticleGuns/interface/FlatRandomPtGunProducer2.h"
#include "EMTFTools/ParticleGuns/interface/FlatRandomLLPGunProducer2.h"
#include "EMTFTools/ParticleGuns/interface/FlatRandomBeamHaloGunProducer2.h"
#include "EMTFTools/ParticleGuns/interface/FlatRandomTau3MuGunProducer2.h"

using edm::FlatRandomPtGunProducer2;
using edm::FlatRandomLLPGunProducer2;
using edm::FlatRandomBeamHaloGunProducer2;
using edm::FlatRandomTau3MuGunProducer2;
DEFINE_FWK_MODULE(FlatRandomPtGunProducer2);
DEFINE_FWK_MODULE(FlatRandomLLPGunProducer2);
DEFINE_FWK_MODULE(FlatRandomBeamHaloGunProducer2);
DEFINE_FWK_MODULE(FlatRandomTau3MuGunProducer2);

