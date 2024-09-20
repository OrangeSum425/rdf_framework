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
                                       std::string name, axis_non_uniform x,
                                       nlohmann::json &cut,
                                       std::string_view wname) {
  auto node = do_cut(noderaw, cut);
  auto &var = x.var;
  auto &bin_edges = x.bin_edges;
  auto nbins = bin_edges.size() - 1;
  auto hist = wname.empty()
                  ? node.Histo1D(ROOT::RDF::TH1DModel{name.data(), var.data(),
                                                      static_cast<int>(nbins),
                                                      bin_edges.data()},
                                 var)
                  : node.Histo1D(ROOT::RDF::TH1DModel{name.data(), var.data(),
                                                      static_cast<int>(nbins),
                                                      bin_edges.data()},
                                 var, wname);
  return hist;
}

ROOT::RDF::RResultPtr<TH2D>
draw_hists_2d(ROOT::RDF::RNode rawnode, std::string name, axis_non_uniform x,
              axis_non_uniform y, nlohmann::json &cut, std::string_view wname) {
  auto node = do_cut(rawnode, cut);
  auto &varx = x.var;
  auto &&xbin_edges = x.bin_edges;
  auto &vary = y.var;
  auto &&ybin_edges = y.bin_edges;
  // auto &ymin = y.min;
  // auto &ymax = y.max;
  // auto &nbinsy = y.bins;
  auto hist =
      wname.empty()
          ? node.Histo2D(
                ROOT::RDF::TH2DModel{
                    name.data(), varx.data(),
                    static_cast<int>(xbin_edges.size() - 1), xbin_edges.data(),
                    static_cast<int>(ybin_edges.size() - 1), ybin_edges.data()},
                varx, vary)

          : node.Histo2D(
                ROOT::RDF::TH2DModel{
                    name.data(), varx.data(),
                    static_cast<int>(xbin_edges.size() - 1), xbin_edges.data(),
                    static_cast<int>(ybin_edges.size() - 1), ybin_edges.data()},
                varx, vary, wname);
  return hist;
}

ROOT::RDF::RResultPtr<TH3D>
draw_hists_3d(ROOT::RDF::RNode rawnode, std::string name, axis_non_uniform x,
              axis_non_uniform y, axis_non_uniform z, nlohmann::json &cut,
              std::string_view wname) {
  auto node = do_cut(rawnode, cut);
  auto &varx = x.var;
  auto &&xbin_edges = x.bin_edges;
  auto &vary = y.var;
  auto &&ybin_edges = y.bin_edges;
  auto &varz = z.var;
  auto &&zbin_edges = z.bin_edges;
  auto hist =
      wname.empty()
          ? node.Histo3D(
                ROOT::RDF::TH3DModel{
                    name.data(), varx.data(),
                    static_cast<int>(xbin_edges.size() - 1), xbin_edges.data(),
                    static_cast<int>(ybin_edges.size() - 1), ybin_edges.data(),
                    static_cast<int>(zbin_edges.size() - 1), zbin_edges.data()},
                varx, vary, varz)

          : node.Histo3D(
                ROOT::RDF::TH3DModel{
                    name.data(), varx.data(),
                    static_cast<int>(xbin_edges.size() - 1), xbin_edges.data(),
                    static_cast<int>(ybin_edges.size() - 1), ybin_edges.data(),
                    static_cast<int>(zbin_edges.size() - 1), zbin_edges.data()},
                varx, vary, varz, wname);
  return hist;
}


axis_non_uniform get_axis(axis &ax) {
  return std::visit(
      [](auto &&arg) -> axis_non_uniform {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, axis_uniform>) {
          axis_uniform &a = arg;
          auto min = a.min;
          auto max = a.max;
          auto bins = a.bins;
          auto bin_width = (max - min) / bins;
          std::vector<double> bin_edges(bins + 1);
          for (int i = 0; i <= bins; i++) {
            bin_edges[i] = min + i * bin_width;
          }
          return axis_non_uniform{.var = a.var, .bin_edges = bin_edges};
        } else if constexpr (std::is_same_v<T, axis_non_uniform>) {
          return arg;
        }
      },
      ax);
}

ROOT::RDF::RResultPtr<TH1> draw_hists_nd(ROOT::RDF::RNode rawnode,
                                         std::string name,
                                         std::vector<axis> axes,
                                         nlohmann::json &cut,
                                         std::string_view wname) {
  switch (axes.size()) {
  case 1:
    return draw_hists(rawnode, name, get_axis(axes[0]), cut, wname);
  case 2:
    return draw_hists_2d(rawnode, name, get_axis(axes[0]), get_axis(axes[1]),
                         cut, wname);
  case 3:
    return draw_hists_3d(rawnode, name, get_axis(axes[0]), get_axis(axes[1]),
                         get_axis(axes[2]), cut, wname);
  default:
    throw std::runtime_error("THnD not implemented yet");
    // the code of drawing the THnD is ready
    // but need to rethink handling of normalization
    // since THnD is not inherited from TH1
    // return draw_hists_nd_imp(rawnode, name, axes, cut, wname);
  }
}