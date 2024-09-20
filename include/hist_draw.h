#pragma once
#include "TH1.h"
#include "TH2.h"
#include <ROOT/RDataFrame.hxx>
#include <nlohmann/json.hpp>
#include <variant>

struct axis_uniform {
  std::string var;
  double min;
  double max;
  int bins;
};

struct axis_non_uniform {
  std::string var;
  std::vector<double> bin_edges;
};

using axis = std::variant<axis_uniform, axis_non_uniform>;

ROOT::RDF::RResultPtr<TH1D> draw_hists(ROOT::RDF::RNode node, std::string name,
                                       axis_non_uniform x, nlohmann::json &cut,
                                       std::string_view wname = "");

ROOT::RDF::RResultPtr<TH2D> draw_hists_2d(ROOT::RDF::RNode node,
                                          std::string name, axis_non_uniform x,
                                          axis_non_uniform y,
                                          nlohmann::json &cut,
                                          std::string_view wname = "");

ROOT::RDF::RResultPtr<TH3D>
draw_hists_3d(ROOT::RDF::RNode node, std::string name, axis_non_uniform x,
              axis_non_uniform y, axis_non_uniform z, nlohmann::json &cut,
              std::string_view wname = "");

ROOT::RDF::RNode do_cut(ROOT::RDF::RNode node, nlohmann::json &cut);

ROOT::RDF::RResultPtr<TH1> draw_hists_nd(ROOT::RDF::RNode rawnode,
                                         std::string name,
                                         std::vector<axis> axes,
                                         nlohmann::json &cut,
                                         std::string_view wname);

axis_non_uniform get_axis(axis &ax);