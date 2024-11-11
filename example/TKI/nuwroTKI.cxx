#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TObjString.h"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <common.h>
#include <tki_general.h>

event::channel getmode_nuwro(TObjString &code) {
  const int mode = code.GetString().Atoi();
  // cout << "mode = " << mode << endl;
  switch (mode) {
  case 1:
    return event::channel::QE;
    break;
  case 11:
    return event::channel::RES;
    break;
  case 26:
    return event::channel::DIS;
    break;
  case 2:
    return event::channel::MEC;
    break;
  case 16:
  default:
    return event::channel::Other;
    break;
  }
}

ROOT::RDF::RNode NuWro_RDF_setup_event(ROOT::RDF::RNode df) {
  return df
      .Define("event",
              [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                 ROOT::RVec<int> &StdHepStatus, ROOT::RVec<double> &StdHepP4_,
                 TObjString &EvtCode) {
                double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                event e{};
                if (getmode_nuwro(EvtCode) == event::channel::Other) {
                  return e;
                }
                e.set_mode(getmode_nuwro(EvtCode));
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
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"})
      .Filter([](event &e) { return e.get_mode() != event::channel::Other; },
              {"event"});
}

class pre : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return MINERvAGFS_general(NuWro_RDF_setup_event(df));
  }
};
 
class normalize_factor_CC : public NormalizeI {
public:
  double operator()(ROOT::RDF::RNode df) override {
    auto mean = df.Mean("EvtWght");
    auto num = df.Count();
    auto norm = mean.GetValue() / num.GetValue();
    norm/= 1e38;
    std::cout << "Cross section per nucleon: " << norm << std::endl;
    return norm;
  }
};

REGISTER_NORMALIZE(normalize_factor_CC);
REGISTER_PROCESS_NODE(pre);

class CCQEPure : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::QE;
        },
        {"event"}, "CUTRES");
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
        {"event"}, "CUTRES");
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
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCMECPure)






class CCQEPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::QE &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCQEPure0pi)

class CCRESPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::RES &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCRESPure0pi)

class CCDISPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::DIS &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCDISPure0pi)

class CCMECPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::MEC &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCMECPure0pi)

class CCotherPure0pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::Other &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211)) == 0;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCotherPure0pi)

class CCRESPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::RES &&
                 (e.count_out(211) + e.count_out(111) + e.count_out(-211))  == 1;
        },
        {"event"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCRESPure1pi)

class CCRESPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::RES &&
                 e.count_out(111)  > 1;
        },
        {"event"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCRESPureMpi)

class CCDISPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::DIS &&
                 e.count_out(111)  == 1;
        },
        {"event"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCDISPure1pi)

class CCDISPureMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::DIS &&
                 e.count_out(111)  > 1;
        },
        {"event"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCDISPureMpi)
/*
class CCRESPYTHIA : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e, int flag_delta) {
          return e.get_mode() == event::channel::RES && !flag_delta;
        },
        {"event", "flag_delta"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCRESPYTHIA)

class CCRESPYTHIA1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e, int flag_delta) {
          return e.get_mode() == event::channel::RES && !flag_delta &&
                 e.count_out(111)  == 1;
        },
        {"event", "flag_delta"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCRESPYTHIA1pi)


class CCRESPYTHIAMpi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e, int flag_delta) {
          return e.get_mode() == event::channel::RES && !flag_delta &&
                 e.count_out(111)  > 1;
        },
        {"event", "flag_delta"}, "CUTRES");
  }
};

REGISTER_PROCESS_NODE(CCRESPYTHIAMpi)*/


class CCQEPure1pi : public ProcessNodeI {
public:
  ROOT::RDF::RNode operator()(ROOT::RDF::RNode df) override {
    return df.Filter(
        [](event &e) {
          return e.get_mode() == event::channel::QE &&
                 e.count_out(111)  == 1;
        },
        {"event"}, "CUTRES");
  }
};
REGISTER_PROCESS_NODE(CCQEPure1pi)