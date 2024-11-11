#include "common.h"
#include <ROOT/RDF/HistoModels.hxx>
#include <TF1.h>
#include <TObjString.h>
#include <TSpline.h>
#include <common.h>
#include <event.h>

class GENIE2Event : public ProcessNodeI {
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df
        .Define(
            "GenEvent",
            [](int StdHepN, const ROOT::RVec<double> &StdHepP4,
               const ROOT::RVec<int> &StdHepPdg,
               const ROOT::RVec<int> &StdHepStatus, const TObjString &EvtCode) {
              event e{};
              auto str = EvtCode.GetString();
              if (str.Contains("QES")) {
                e.set_interaction(interaction_channel::qel);
              } else if (str.Contains("RES")) {
                e.set_interaction(interaction_channel::res);
              } else if (str.Contains("DIS")) {
                e.set_interaction(interaction_channel::dis);
              } else if (str.Contains("MEC")) {
                e.set_interaction(interaction_channel::mec);
              } else if (str.Contains("COH")) {
                e.set_interaction(interaction_channel::coh);
              } else {
                e.set_interaction(interaction_channel::unknown);
              }

              e.set_is_cc(str.Contains("CC"));

              for (int i = 0; i < StdHepN; i++) {
                switch (StdHepStatus[i]) {
                case 0:
                  e.add_particle_in(
                      StdHepPdg[i],
                      TLorentzVector{StdHepP4[4 * i + 0], StdHepP4[4 * i + 1],
                                     StdHepP4[4 * i + 2], StdHepP4[4 * i + 3]});
                  break;
                case 1:
                  e.add_particle_out(
                      StdHepPdg[i],
                      TLorentzVector{StdHepP4[4 * i + 0], StdHepP4[4 * i + 1],
                                     StdHepP4[4 * i + 2], StdHepP4[4 * i + 3]});
                  break;
                case 11:
                  e.set_hit_nucleon(
                      TLorentzVector{StdHepP4[4 * i + 0], StdHepP4[4 * i + 1],
                                     StdHepP4[4 * i + 2], StdHepP4[4 * i + 3]});
                default:
                  break;
                }
              }

              return e;
            },
            {"StdHepN", "StdHepP4", "StdHepPdg", "StdHepStatus", "EvtCode"})
        // .Filter(
        //     [](const event &GenEvent) {
        //       return GenEvent.get_interaction() !=
        //                  interaction_channel::unknown &&
        //              GenEvent.get_interaction() != interaction_channel::coh;
        //     },
        //     {"GenEvent"})
        .Define("cc_flag", [](const event &e) { return e.get_is_cc(); },
                {"GenEvent"});
  };

  void configure(const nlohmann::json &) override {};
};

double get_fuxint(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  // spline->SaveAs("wrong.root");
  TSpline3 sp("sp", spline);
  TF1 func(
      "spline", [&](double *x, double *) { return sp.Eval(*x); }, 0,
      h_rate->GetXaxis()->GetXmax(), 0);
  for (int ii = 1; ii <= h_rate->GetNbinsX(); ii++) {
    double bin_c = h_rate->GetBinContent(ii);
    double bin_up = h_rate->GetXaxis()->GetBinUpEdge(ii);
    double bin_low = h_rate->GetXaxis()->GetBinLowEdge(ii);
    double bin_width = bin_up - bin_low;
    if (bin_c < 1 || func.Integral(bin_low, bin_up) == 0) {
      continue;
    }
    fluxint += bin_c / func.Integral(bin_low, bin_up) * bin_width;
  }
  return fluxint;
}

constexpr double average_A() {
  // mass fraction
  const double kC = 0.8786909745361763;
  const double kH = 0.12130902546382372;
  // mass number
  const double mC = 12.;
  const double mH = 1.007825;
  return (12 * kC * mH + kH * mC) / (kC * mH + kH * mC);
}

class MIX_NORMALIZE : public NormalizeI {
  // this will only work for CH
public:
  double operator()(ROOT::RDF::RNode df) override {
    TFile root_file{filename.c_str(), "READ"};
    auto nu_mu_C12 =
        static_cast<TGraph *>(root_file.Get("nu_mu_C12/tot_cc")->Clone());
    auto nu_mu_bar_C12 =
        static_cast<TGraph *>(root_file.Get("nu_mu_bar_C12/tot_cc")->Clone());
    auto nu_mu_H1 =
        static_cast<TGraph *>(root_file.Get("nu_mu_H1/tot_cc")->Clone());
    auto nu_mu_bar_H1 =
        static_cast<TGraph *>(root_file.Get("nu_mu_bar_H1/tot_cc")->Clone());

    auto nu_e_C12 =
        static_cast<TGraph *>(root_file.Get("nu_e_C12/tot_cc")->Clone());
    auto nu_e_bar_C12 =
        static_cast<TGraph *>(root_file.Get("nu_e_bar_C12/tot_cc")->Clone());
    auto nu_e_H1 =
        static_cast<TGraph *>(root_file.Get("nu_e_H1/tot_cc")->Clone());
    auto nu_e_bar_H1 =
        static_cast<TGraph *>(root_file.Get("nu_e_bar_H1/tot_cc")->Clone());

    ROOT::RDF::TH1DModel h_model{"", "", 20, 0, 20.};
    auto dfcc = df.Filter(
        [](const TObjString &EvtCode) {
          return EvtCode.GetString().Contains("CC");
        },
        {"EvtCode"});
    auto get_hist_neutrinoE_cc = [&](int neutrino, int nucleus) {
      return dfcc
          .Filter(
              [=](const ROOT::RVec<int> &StdHepPdg) {
                return StdHepPdg[0] == neutrino && StdHepPdg[1] == nucleus;
              },
              {"StdHepPdg"})
          .Define(
              "neutrinoE",
              [](const ROOT::RVec<double> &StdHepP4) { return StdHepP4[3]; },
              {"StdHepP4"})
          .Histo1D(h_model, "neutrinoE");
    };    

    auto h_nu_mu_C12 = get_hist_neutrinoE_cc(14, 1000060120);
    auto h_nu_mu_bar_C12 = get_hist_neutrinoE_cc(-14, 1000060120);
    auto h_nu_mu_H1 = get_hist_neutrinoE_cc(14, 2212);
    auto h_nu_mu_bar_H1 = get_hist_neutrinoE_cc(-14, 2212);

    auto h_nu_e_C12 = get_hist_neutrinoE_cc(12, 1000060120);
    auto h_nu_e_bar_C12 = get_hist_neutrinoE_cc(-12, 1000060120);
    auto h_nu_e_H1 = get_hist_neutrinoE_cc(12, 2212);
    auto h_nu_e_bar_H1 = get_hist_neutrinoE_cc(-12, 2212);

    auto event_count = dfcc.Count().GetValue();
    auto total_fluxint =
        get_fuxint(h_nu_mu_C12.GetPtr(), nu_mu_C12)* 12+
        get_fuxint(h_nu_mu_bar_C12.GetPtr(), nu_mu_bar_C12) * 12 +
        get_fuxint(h_nu_mu_H1.GetPtr(), nu_mu_H1) * 1.6 +
        get_fuxint(h_nu_mu_bar_H1.GetPtr(), nu_mu_bar_H1) * 1.6 +
        get_fuxint(h_nu_e_C12.GetPtr(), nu_e_C12) * 12 +
        get_fuxint(h_nu_e_bar_C12.GetPtr(), nu_e_bar_C12) * 12+
        get_fuxint(h_nu_e_H1.GetPtr(), nu_e_H1) * 1.6 +
        get_fuxint(h_nu_e_bar_H1.GetPtr(), nu_e_bar_H1) * 1.6;
        
    auto xsec_per_nucleon = event_count / total_fluxint;
    // per nucleus to per nucleon
    auto res = xsec_per_nucleon;
    std::cout << "Total event count: " << event_count << std::endl;
    std::cout << "Total fluxint: " << total_fluxint << std::endl;
    std::cout << "Cross section per nucleon: " << xsec_per_nucleon << std::endl;
    return res / event_count;
    //return res;
  }

  void configure(const nlohmann::json &conf) override {
    filename = conf["filename"];
  }

private:
  std::string filename;
};

REGISTER_PROCESS_NODE(GENIE2Event)
REGISTER_NORMALIZE(MIX_NORMALIZE)