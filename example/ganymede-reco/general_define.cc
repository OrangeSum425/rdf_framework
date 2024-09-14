#include "TLorentzVector.h"
#include <common.h>
#include <event.h>

constexpr double m_N = (0.93827210 + 0.93956540) / 2;

class ANALYSIS : public ProcessNodeI {
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df
        .Define("Q2",
                [](const event &GenEvent) {
                  auto lvq =
                      GenEvent.get_neutrino() - GenEvent.get_primary_lepton();
                  return -lvq.M2();
                },
                {"GenEvent"})
        .Define("W",
                [](const event &GenEvent) {
                  auto p4_hadronic_system = GenEvent.get_nucleon() +
                                            GenEvent.get_neutrino() -
                                            GenEvent.get_primary_lepton();
                  return p4_hadronic_system.M();
                },
                {"GenEvent"})
        .Define("xBj",
                [](const event &GenEvent, double Q2) {
                  auto lvq =
                      GenEvent.get_neutrino() - GenEvent.get_primary_lepton();
                  auto initN = GenEvent.get_nucleon();
                  return Q2 / (2 * initN.Dot(lvq));
                },
                {"GenEvent", "Q2"})
        .Define("W_rest",
                [](const event &GenEvent) {
                  auto p4_hadronic_system = TLorentzVector{0, 0, 0, m_N} +
                                            GenEvent.get_neutrino() -
                                            GenEvent.get_primary_lepton();
                  return p4_hadronic_system.M();
                },
                {"GenEvent"})
        .Define("xBj_rest",
                [](const event &GenEvent, double Q2) {
                  auto lvq =
                      GenEvent.get_neutrino() - GenEvent.get_primary_lepton();
                  auto initN = TLorentzVector{0, 0, 0, m_N};
                  return Q2 / (2 * initN.Dot(lvq));
                },
                {"GenEvent", "Q2"})
        .Define("y",
                [](const event &GenEvent) {
                  auto enu = GenEvent.get_neutrino().E();
                  auto el = GenEvent.get_primary_lepton().E();
                  return 1 - el / enu;
                },
                {"GenEvent"})
        .Define(
            "lepton_p4",
            [](const event &GenEvent) { return GenEvent.get_primary_lepton(); },
            {"GenEvent"})
        .Define(
            "Enu",
            [](const event &GenEvent) { return GenEvent.get_neutrino().E(); },
            {"GenEvent"})
        .Define("El", [](const TLorentzVector &p4) { return p4.E(); },
                {"lepton_p4"})
        .Define("q0", [](double enu, double el) { return enu - el; },
                {"Enu", "El"})
        .Define("neutrino_p4",
                [](const event &GenEvent) { return GenEvent.get_neutrino(); },
                {"GenEvent"})
        .Define("npi0",
                [](const event &GenEvent) {
                  return GenEvent.count_particle_out(111);
                },
                {"GenEvent"})
        .Define("npi_charged",
                [](const event &GenEvent) {
                  return GenEvent.count_particle_out(211) +
                         GenEvent.count_particle_out(-211);
                },
                {"GenEvent"})
        .Define("n_neutron",
                [](const event &GenEvent) {
                  return GenEvent.count_particle_out(2112);
                },
                {"GenEvent"})
        .Define("n_neutron_d",
                [](size_t n_neutron) { return static_cast<double>(n_neutron); },
                {"n_neutron"}) // just for the sake of histogramming
                               // (RDF::Histo1D wants double)
        .Define("is_qel",
                [](const event &GenEvent) {
                  return GenEvent.get_interaction() == interaction_channel::qel;
                },
                {"GenEvent"});
  }

  void configure(const nlohmann::json &) override {};
};

REGISTER_PROCESS_NODE(ANALYSIS)
