#ifndef FlatRandomPtGunProducer2_H
#define FlatRandomPtGunProducer2_H

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

    class FlatRandomPtGunProducer2 : public BaseFlatGunProducer2 {
        public:
            explicit FlatRandomPtGunProducer2(const ParameterSet&);
            ~FlatRandomPtGunProducer2() override;

            void produce(Event&, const EventSetup&) override;

        private:
            double min_pt_;
            double max_pt_;
            double min_dxy_;
            double max_dxy_;
            double min_one_over_pt_;
            double max_one_over_pt_;
            double min_one_over_dxy_;
            double max_one_over_dxy_;
            double max_lxy_;
            bool random_charge_en_;
            std::string pt_spectrum_;
            std::string vertex_spectrum_;
            
            // Utils
            double calcEta(
                    const double&, 
                    const double&,
                    const double&,
                    const double&,
                    const double&,
                    const double&
            ) const;

            double calcEtaSt2(
                    const double&, 
                    const double&,
                    const double&,
                    const double&,
                    const double&,
                    const double&
            ) const;

            HepMC::FourVector newVertexWithD0(
                    CLHEP::HepRandomEngine*, 
                    const double&, 
                    const double&,
                    const double&
            ) const;
    };

}  // namespace edm

#endif
