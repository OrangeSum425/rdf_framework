#include "TF1.h"
#include "TRandom3.h"
#include "common.h"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TSpline.h>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <tki_general.h>
#include <type_traits>
#include <vector>

event::channel get_mode_genie(const TObjString &code) {
  if (code.GetString().Contains("QES")) {
    return event::channel::QE;
  } else if (code.GetString().Contains("RES")) {
    return event::channel::RES;
  } else if (code.GetString().Contains("DIS")) {
    return event::channel::DIS;
  } else if (code.GetString().Contains("MEC")) {
    return event::channel::MEC;
  }
  return event::channel::Other;
}

ROOT::RDF::RNode GENIE_RDF_setup_event(ROOT::RDF::RNode df) {
  return df
      .Filter(
          [](TObjString &EvtCode) {
            auto str = EvtCode.GetString();
            return str.Contains("CC");
          },
          {"EvtCode"})
      .Define("neutrinoE",
              [](int StdHepN, ROOT::RVec<double> &StdHepP4,
                 ROOT::RVec<int> &StdHepPdg, ROOT::RVec<int> &StdHepStatus) {
                for (int i = 0; i < StdHepN; i++) {
                  if ((abs(StdHepPdg[i]) == 12 || abs(StdHepPdg[i]) == 14 ||
                       abs(StdHepPdg[i]) == 16) &&
                      StdHepStatus[i] == 0) {
                    return StdHepP4[4 * i + 3];
                  }
                }
                return 0.;
              },
              {"StdHepN", "StdHepP4", "StdHepPdg", "StdHepStatus"})
      .Define("event",
              [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                 ROOT::RVec<int> &StdHepStatus, ROOT::RVec<double> &StdHepP4_,
                 TObjString &EvtCode) {
                double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                event e{};
                e.set_mode(get_mode_genie(EvtCode));
                for (int i = 0; i < StdHepN; ++i) {
                  auto pdg = StdHepPdg[i];
                  if (StdHepPdg[i] == 1000000010) {
                    pdg = 2112;
                  }
                  switch (StdHepStatus[i]) {
                  case 0:
                    e.add_particle_in(
                        pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                  case 1:
                    e.add_particle_out(
                        pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                  case 2: {
                    e.add_particle_nofsi(
                        pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
                  } break;
                  default:
                    break;
                  }
                }
                return e;
              },
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"});
}





class pre : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return MINERvAGFS_general(GENIE_RDF_setup_event(df));
  }
};


double get_fuxint(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  TSpline3 sp("sp", spline);
  TF1 func("spline", [&](double *x, double *) { return sp.Eval(*x); }, 0,
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

// Shared class for normalization across nodes
class normalize_factor_CC : public NormalizeI {
public:
  double operator()(ROOT::RDF::RNode df) override {
    TFile root_file{filename.c_str(), "READ"};
    auto nu_mu_C12 = static_cast<TGraph *>(root_file.Get("nu_mu_C12/tot_cc")->Clone());

    ROOT::RDF::TH1DModel h_model{"", "", 20, 0, 20};
    auto dfcc = df.Filter([](const TObjString &EvtCode) {
      return EvtCode.GetString().Contains("CC");
    }, {"EvtCode"});
    
    auto get_hist_neutrinoE_cc = [&](int neutrino, int nucleus) {
      return dfcc.Filter([=](const ROOT::RVec<int> &StdHepPdg) {
        return StdHepPdg[0] == neutrino && StdHepPdg[1] == nucleus;
      }, {"StdHepPdg"})
      .Define("neutrinoE", [](const ROOT::RVec<double> &StdHepP4) {
        return StdHepP4[3];
      }, {"StdHepP4"}).Histo1D(h_model, "neutrinoE");
    };

    auto h_nu_mu_C12 = get_hist_neutrinoE_cc(14, 1000060120);
    auto event_count = dfcc.Count().GetValue();
    auto total_fluxint = get_fuxint(h_nu_mu_C12.GetPtr(), nu_mu_C12) * 13;
    
    auto xsec_per_nucleon = event_count / total_fluxint;
    xsec_per_nucleon /= 1e38;

    std::cout << "Total event count: " << event_count << std::endl;
    std::cout << "Total fluxint: " << total_fluxint << std::endl;
    std::cout << "Cross section per nucleon (divided by 10^38): " << xsec_per_nucleon << std::endl;

    return xsec_per_nucleon / event_count;
  }

  void configure(const nlohmann::json &conf) override {
    filename = conf["filename"];
  }

private:
  std::string filename;
};

REGISTER_NORMALIZE(normalize_factor_CC);
REGISTER_PROCESS_NODE(pre);

// ProcessNode classes for filtering modes
class CCQEPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::QE &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCQEPure0pi)

class CCQEPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::QE &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCQEPure1pi)

class CCRESPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCRESPure0pi)

class CCDISPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCDISPure0pi)

class CCMECPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::MEC &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCMECPure0pi)

class CCRESPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCRESPure1pi)

class CCDISPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCDISPure1pi)

class CCRESPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) > 1;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCRESPureMpi)

class CCDISPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
             (e.count_out(211) + e.count_out(111) + e.count_out(-211)) > 1;
    }, {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCDISPureMpi)
