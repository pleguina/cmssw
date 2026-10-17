#ifndef IOMC_BaseEvtVtxGenerator2_H
#define IOMC_BaseEvtVtxGenerator2_H

#include <string>
#include <memory>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "TMatrixD.h"

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

// Class
class BaseEvtVtxGenerator2 : public edm::stream::EDProducer<> {
    public:
        // ctor & dtor
        explicit BaseEvtVtxGenerator2(const edm::ParameterSet&);
        ~BaseEvtVtxGenerator2() override;

        void produce(edm::Event&, const edm::EventSetup&) override;

        virtual HepMC::FourVector newVertex(CLHEP::HepRandomEngine*) const = 0;

        virtual TMatrixD const* GetInvLorentzBoost() const = 0;

    protected:
        edm::EDGetTokenT<edm::HepMCProduct> source_token_;
};

#endif
