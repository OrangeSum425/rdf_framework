#include "hist_draw.h"
#include "get_callable.h"
#include <THn.h>
#include <TObject.h>
#include <string>
#include <string_view>
#include <variant>

ROOT::RDF::RNode do_cut(ROOT::RDF::RNode node, nlohmann::json &cut) {
  auto cut_func =
      cut.empty() ? std::make_unique<noop>() : get_node_process_callable(cut);
  if (cut_func)
    return (*cut_func)(node);
  else {
    auto cutstring = cut.get<std::string>();
    // std::cout << "Failed to parse cut: " << cutstring << " as process node"
    //           << " passing it as string cut instead" << std::endl;
    return node.Filter(cut.get<std::string>());
  }
}

ROOT::RDF::RResultPtr<TH1D> draw_hists(ROOT::RDF::RNode noderaw,
                                       std::string name, axis x,
                                       nlohmann::json &cut,
                                       std::string_view wname) {
  auto node = do_cut(noderaw, cut);
  auto &var = x.var;
  auto &xmin = x.min;
  auto &xmax = x.max;
  auto &nbins = x.bins;
  auto hist = wname.empty()
                  ? node.Histo1D(ROOT::RDF::TH1DModel{name.data(), var.data(),
                                                      nbins, xmin, xmax},
                                 var)
                  : node.Histo1D(ROOT::RDF::TH1DModel{name.data(), var.data(),
                                                      nbins, xmin, xmax},
                                 var, wname);
  return hist;
}

ROOT::RDF::RResultPtr<TH2D> draw_hists_2d(ROOT::RDF::RNode rawnode,
                                          std::string name, axis x, axis y,
                                          nlohmann::json &cut,
                                          std::string_view wname) {
  auto node = do_cut(rawnode, cut);
  auto &varx = x.var;
  auto &xmin = x.min;
  auto &xmax = x.max;
  auto &nbinsx = x.bins;
  auto &vary = y.var;
  auto &ymin = y.min;
  auto &ymax = y.max;
  auto &nbinsy = y.bins;
  auto hist =
      wname.empty()
          ? node.Histo2D(ROOT::RDF::TH2DModel{name.data(), varx.data(), nbinsx,
                                              xmin, xmax, nbinsy, ymin, ymax},
                         varx, vary)

          : node.Histo2D(ROOT::RDF::TH2DModel{name.data(), varx.data(), nbinsx,
                                              xmin, xmax, nbinsy, ymin, ymax},
                         varx, vary, wname);
  return hist;
}

ROOT::RDF::RResultPtr<TH3D> draw_hists_3d(ROOT::RDF::RNode rawnode,
                                          std::string name, axis x, axis y,
                                          axis z, nlohmann::json &cut,
                                          std::string_view wname) {
  auto node = do_cut(rawnode, cut);
  auto &varx = x.var;
  auto &xmin = x.min;
  auto &xmax = x.max;
  auto &nbinsx = x.bins;
  auto &vary = y.var;
  auto &ymin = y.min;
  auto &ymax = y.max;
  auto &nbinsy = y.bins;
  auto &varz = z.var;
  auto &zmin = z.min;
  auto &zmax = z.max;
  auto &nbinsz = z.bins;
  auto hist =
      wname.empty()
          ? node.Histo3D(ROOT::RDF::TH3DModel{name.data(), varx.data(), nbinsx,
                                              xmin, xmax, nbinsy, ymin, ymax,
                                              nbinsz, zmin, zmax},
                         varx, vary, varz)

          : node.Histo3D(ROOT::RDF::TH3DModel{name.data(), varx.data(), nbinsx,
                                              xmin, xmax, nbinsy, ymin, ymax,
                                              nbinsz, zmin, zmax},
                         varx, vary, varz, wname);
  return hist;
}

ROOT::RDF::RResultPtr<THn> draw_hists_nd_imp(ROOT::RDF::RNode rawnode,
                                             std::string name,
                                             std::vector<axis> axes,
                                             nlohmann::json &cut,
                                             std::string_view wname) {
  auto node = do_cut(rawnode, cut);
  auto dim = axes.size();
  std::vector<double> bin_low, bin_up;
  std::vector<int> bins;
  std::vector<std::string> vars;
  for (auto &ax : axes) {
    vars.push_back(ax.var);
    bin_low.push_back(ax.min);
    bin_up.push_back(ax.max);
    bins.push_back(ax.bins);
  }
  if (!wname.empty()) {
    vars.push_back(wname.data());
  }
  return node.HistoND(ROOT::RDF::THnDModel{name.data(), name.data(), (int)dim,
                                           bins, bin_low, bin_up},
                      vars);
}

// std::variant<ROOT::RDF::RResultPtr<TH1>, ROOT::RDF::RResultPtr<THn>>
ROOT::RDF::RResultPtr<TH1> draw_hists_nd(ROOT::RDF::RNode rawnode,
                                         std::string name,
                                         std::vector<axis> axes,
                                         nlohmann::json &cut,
                                         std::string_view wname) {
  switch (axes.size()) {
  case 1:
    return draw_hists(rawnode, name, axes[0], cut, wname);
  case 2:
    return draw_hists_2d(rawnode, name, axes[0], axes[1], cut, wname);
  case 3:
    return draw_hists_3d(rawnode, name, axes[0], axes[1], axes[2], cut, wname);
  default:
    throw std::runtime_error("THnD not implemented yet");
    // the code of drawing the THnD is ready
    // but need to rethink handling of normalization
    // since THnD is not inherited from TH1
    // return draw_hists_nd_imp(rawnode, name, axes, cut, wname);
  }
}