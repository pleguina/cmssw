#ifndef FlatRandomTau3MuGunProducer2_H
#define FlatRandomTau3MuGunProducer2_H

#include "Math/LorentzVector.h"

#include "EMTFTools/ParticleGuns/interface/BaseFlatGunProducer2.h"

// Forward Declare
namespace HepMC {
    class FourVector;
}

namespace CLHEP {
    class HepRandomEngine;
}

namespace edm {
    class HepMCProduct;
}

namespace edm {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVectorD;

    class FlatRandomTau3MuGunProducer2 : public BaseFlatGunProducer2 {

        public:
            explicit FlatRandomTau3MuGunProducer2(const ParameterSet&);
            ~FlatRandomTau3MuGunProducer2() override;

            void produce(Event&, const EventSetup&) override;

        private:
            double min_pt_Ds_;
            double max_pt_Ds_;
            double min_invpt_Ds_;
            double max_invpt_Ds_;
            bool random_charge_en_;

            void shoot_tau(
                    CLHEP::HepRandomEngine*, 
                    const double&, const double&, const double&, 
                    XYZTLorentzVectorD&, XYZTLorentzVectorD&, XYZTLorentzVectorD&
            ) const;

            void decay_meson(
                    CLHEP::HepRandomEngine*, 
                    const double&, const double&, const double&, 
                    const XYZTLorentzVectorD&, XYZTLorentzVectorD&, 
                    XYZTLorentzVectorD&, XYZTLorentzVectorD&
            ) const;

            void decay_tau(
                    CLHEP::HepRandomEngine*, 
                    const double&, const double&, const double&, 
                    const XYZTLorentzVectorD&, XYZTLorentzVectorD&, 
                    XYZTLorentzVectorD&, XYZTLorentzVectorD&, XYZTLorentzVectorD&
            ) const;
    };

}  // namespace edm

#endif
