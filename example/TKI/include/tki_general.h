#pragma once

#include <TLorentzVector.h>
#include <map>

#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>

template <typename T> class eq_range {
private:
  T pair;

public:
  template <typename T_> eq_range(T_ &&p) : pair(p) {}
  decltype(auto) begin() { return pair.first; };
  decltype(auto) end() { return pair.second; };
};
template <typename T_> eq_range(T_ &&) -> eq_range<std::remove_reference_t<T_>>;
class event {
public:
  enum class channel { QE, RES, DIS, MEC, Other };
  bool passesCC0piSelection() const;
  bool passesCCpi0Selection() const;

private:
  std::unordered_multimap<int, TLorentzVector> in_particles{};
  std::unordered_multimap<int, TLorentzVector> out_particles{};
  // std::set<int> pdg_list_out{};
  std::map<int, size_t> pdg_list_out{};
  std::unordered_multimap<int, TLorentzVector> nofsi_particles{};
  double weight{1};
  channel mode;
  TLorentzVector primary_proton{}, spectator_proton{};
  bool found_primary_proton{};
  TLorentzVector primary_lepton_out{};
  TLorentzVector primary_lepton_in{};
  int leptoncharge{};
  bool found_lepton_out{false};

public:
  event() {}
  ~event();
  void setprimaryP(const TLorentzVector &p) { primary_proton = p; }
  const TLorentzVector &getprimaryP() const { return primary_proton; }
  const TLorentzVector &getspectatorP() const { return spectator_proton; }
  double getQ2() const;
  double getW() const;
  void add_particle_in(int id, const TLorentzVector &p4);
  void add_particle_out(int id, const TLorentzVector &p4);
  void add_particle_nofsi(int id, const TLorentzVector &p4);
  std::pair<double, double> get_q2_w() const;
  double get_enu() const;
  TLorentzVector get_leading_proton() const;
  TLorentzVector get_leading_out(int) const;
  TLorentzVector get_leading_nofsi(int) const;
  auto get_particle_out(int pdgid) {
    return eq_range{out_particles.equal_range(pdgid)};
  }
  auto get_particle_in(int pdgid) {
    return eq_range{in_particles.equal_range(pdgid)};
  }
  auto get_particle_nofsi(int pdgid) {
    return eq_range{nofsi_particles.equal_range(pdgid)};
  }
  size_t count_particle_out(int pdgid) const noexcept;
  size_t count_particle_nofsi(int pdgid) const noexcept;
  void set_mode(channel);
  channel get_mode() const;
  void set_weight(double);
  double get_weight() const;
  const std::unordered_multimap<int, TLorentzVector> &get_particle_out() const {
    return out_particles;
  };
  const std::unordered_multimap<int, TLorentzVector> &get_particle_in() const {
    return in_particles;
  };
  const std::unordered_multimap<int, TLorentzVector> &
  get_particle_nofsi() const {
    return nofsi_particles;
  };
  size_t count_out(int) const noexcept;
  void set_nod(const std::array<int, 18> &nod);
  bool is_no_pion_interaction() const;
  bool is_true_elastic() const;
  const std::array<int, 18> &get_nod() const;
  int get_pion_interaction_count() const;
  double getW_nofsi() const;
  // void setPrimaryLepton(const TLorentzVector &p) { primary_lepton_out = p; }
  const TLorentzVector &get_primary_lepton_out() { return primary_lepton_out; }
  const TLorentzVector &get_primary_lepton_in() { return primary_lepton_in; }
  double W_rest() const;
  TLorentzVector get_lvq() const;
  std::map<int, size_t> get_pdg_list_out() const { return pdg_list_out; }
};

template <typename T> inline TLorentzVector get_sum(T &&cont, int pdgid) {
  TLorentzVector sum{};
  for (const auto &[id, p4] : eq_range{cont.equal_range(pdgid)}) {
    sum += p4;
  }
  return sum;
}

template <typename T> inline TLorentzVector get_leading(T &&cont, int pdgid) {
  const TLorentzVector *max{nullptr};
  double max_E{-INFINITY};
  for (const auto &[id, p4] : eq_range{cont.equal_range(pdgid)}) {
    if (p4.E() > max_E) {
      max_E = p4.E();
      max = &p4;
    }
  }
  if (max) {
    return *max;
  } else {
    return TLorentzVector{};
  }
}

std::array<double, 8> getCommonTKI(const int targetA, const int targetZ,
                                   const TLorentzVector *tmp4pBeam,
                                   const TLorentzVector *tmp4pScatter,
                                   const TLorentzVector *tmp4pRecoil);

ROOT::RDF::RNode MINERvAGFSPIZERO_f(ROOT::RDF::RNode df);
ROOT::RDF::RNode MINERvAGFS0PI_f(ROOT::RDF::RNode df);

REGISTER_PROCESS_NODE(MINERvAGFS0PI)

ROOT::RDF::RNode MINERvAGFS_general(ROOT::RDF::RNode df);
ROOT::RDF::RNode MINERvAGFS_do_TKI(ROOT::RDF::RNode df);

