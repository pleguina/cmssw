#ifndef FlatRandomLLPGunProducer2_H
#define FlatRandomLLPGunProducer2_H

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

    class FlatRandomLLPGunProducer2 : public BaseFlatGunProducer2 {

        public:
            explicit FlatRandomLLPGunProducer2(const ParameterSet&);
            ~FlatRandomLLPGunProducer2() override;

            void produce(Event&, const EventSetup&) override;

        private:
            double min_m_h_;
            double max_m_h_;
            double min_pt_h_;
            double max_pt_h_;
            double min_invpt_h_;
            double max_invpt_h_;
            double min_ctau_llp_;
            double max_ctau_llp_;
            std::string llp_mass_spectrum_;

            void shoot_llp(
                    CLHEP::HepRandomEngine*, 
                    const double&, const double&, XYZTLorentzVectorD&, 
                    XYZTLorentzVectorD&, XYZTLorentzVectorD&
            ) const;

            void decay_particle(
                    CLHEP::HepRandomEngine*, 
                    const double&, const double&, const double&, 
                    const XYZTLorentzVectorD&, XYZTLorentzVectorD&, 
                    XYZTLorentzVectorD&, XYZTLorentzVectorD&
            ) const;
    };

}  // namespace edm

#endif
