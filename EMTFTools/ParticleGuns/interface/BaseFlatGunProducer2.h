#ifndef BaseFlatGunProducer2_H
#define BaseFlatGunProducer2_H

#include <string>
#include <memory>
#include <vector>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTESSource/interface/HepPDTESSource.h"

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

namespace edm {

    class BaseFlatGunProducer2 : public one::EDProducer<one::WatchRuns, EndRunProducer> {
        public:
            explicit BaseFlatGunProducer2(const ParameterSet&);
            ~BaseFlatGunProducer2() override;

            void beginRun(const edm::Run&, const edm::EventSetup&) override;
            void endRun(const edm::Run&, const edm::EventSetup&) override;
            void endRunProduce(edm::Run&, const edm::EventSetup&) override;

        protected:
            // gun particle(s) characteristics
            std::vector<int> particle_id_list_;
            double min_eta_;
            double max_eta_;
            double min_phi_;
            double max_phi_;

            // the event format itself
            HepMC::GenEvent* event_;

            // HepMC/HepPDT related things
            int verbosity_;
            bool add_anti_particle_en_;

            ESGetToken<HepPDT::ParticleDataTable, PDTRecord> pdg_table_es_token_;
    };

}  // namespace edm

#endif
