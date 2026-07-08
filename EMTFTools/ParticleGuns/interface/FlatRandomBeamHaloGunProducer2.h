#ifndef FlatRandomBeamHaloGunProducer2_H
#define FlatRandomBeamHaloGunProducer2_H

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

    class FlatRandomBeamHaloGunProducer2 : public BaseFlatGunProducer2 {
        public:
            explicit FlatRandomBeamHaloGunProducer2(const ParameterSet&);
            ~FlatRandomBeamHaloGunProducer2() override;

            void produce(Event&, const EventSetup&) override;

        private:
            double min_p_;
            double max_p_;
            double min_vr_;
            double max_vr_;
            double min_vz_;
            double max_vz_;
            double min_tz_;
            double max_tz_;
            bool random_charge_en_;
            std::string p_spectrum_;
    };

}  // namespace edm

#endif
