#include "TF1.h"
#include "TRandom3.h"
#include "common.h"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
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
//#include "common.h"

class CC_selection : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df
        .Filter([](TObjString &str) { return str.GetString().Contains("CC"); },
                {"EvtCode"}, "CC_cut")
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
                {"StdHepN", "StdHepP4", "StdHepPdg", "StdHepStatus"});
  }
};

REGISTER_PROCESS_NODE(CC_selection)


template <typename T>
std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
  TFile root_file{file_path.c_str(), "READ"};
  auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
  assert(objptr);
  return std::unique_ptr<T>{objptr};
}

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


std::pair<double, double> get_xsec_tki(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  // spline->SaveAs("wrong.root");
  TSpline3 sp("sp", spline);
  TF1 func(
      "spline", [&](double *x, double *) { return sp.Eval(*x) / 1e38; }, 0,
      h_rate->GetXaxis()->GetXmax(), 0);
  double integral_total = func.Integral(0, h_rate->GetXaxis()->GetXmax());
  std::cerr << "Total integral : " << integral_total << std::endl;
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
  double event_rate = h_rate->Integral();
  return {event_rate, event_rate / fluxint};
}

class normalize_factor_cc : public NormalizeI {
public:
  double operator()(ROOT::RDF::RNode df) override {
    auto h = (CC_selection{})(df).Histo1D({"", "", 200, 0, 100}, "neutrinoE");
    // auto filename = parameters[0];
    // auto obj_path = parameters[1];
    // auto Z = std::stoi(parameters[2]);
    auto spline_obj = get_object<TGraph>(filename, obj_path);
    auto [tot, xsec] = get_xsec_tki(h.GetPtr(), spline_obj.get());
    std::cerr << "tot: " << tot << " xsec: " << xsec << std::endl;
    xsec *= 1 / ((double)Z);
    return xsec / tot;
  }

  void configure(const nlohmann::json &conf) override {
    filename = conf["filename"];
    obj_path = conf["obj_path"];
    Z = conf["Z"];
  }

private:
  std::string filename;
  std::string obj_path;
  int Z;
};
REGISTER_NORMALIZE(normalize_factor_cc);
REGISTER_PROCESS_NODE(pre);


class CCQEPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::QE &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTQE0PI");
  }
};
REGISTER_PROCESS_NODE(CCQEPure0pi)

class CCQEPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::QE &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTQE1PI");
  }
};
REGISTER_PROCESS_NODE(CCQEPure1pi)

class CCRESPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTRES0PI");
  }
};
REGISTER_PROCESS_NODE(CCRESPure0pi)

class CCDISPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTDIS0PI");
  }
};
REGISTER_PROCESS_NODE(CCDISPure0pi)

class CCMECPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::MEC &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
    }, {"event"}, "CUTMEC0PI");
  }
};
REGISTER_PROCESS_NODE(CCMECPure0pi)

class CCRESPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTRES1PI");
  }
};
REGISTER_PROCESS_NODE(CCRESPure1pi)

class CCDISPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 1;
    }, {"event"}, "CUTDIS1PI");
  }
};
REGISTER_PROCESS_NODE(CCDISPure1pi)

class CCRESPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::RES &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) > 1;
    }, {"event"}, "CUTRESMPI");
  }
};
REGISTER_PROCESS_NODE(CCRESPureMpi)

class CCDISPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter([](event &e) {
      return e.get_mode() == event::channel::DIS &&
              (e.count_out(211) + e.count_out(111) + e.count_out(-211)) > 1;
    }, {"event"}, "CUTDISMPI");
  }
};
REGISTER_PROCESS_NODE(CCDISPureMpi)

class CCQEPure : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::QE;
        },
        {"event"}, "CUTQE");
  }
};
REGISTER_PROCESS_NODE(CCQEPure)

class CCRESPure : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::RES;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCRESPure)

class CCDISPure : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::DIS;
        },
        {"event"}, "CUTDIS");
  }
};
REGISTER_PROCESS_NODE(CCDISPure)

class CCMECPure : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::MEC;
        },
        {"event"}, "CUTMEC");
  }
};
REGISTER_PROCESS_NODE(CCMECPure)
