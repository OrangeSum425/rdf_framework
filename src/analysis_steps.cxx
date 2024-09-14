#include "analysis_steps.h"
#include "THStack.h"
#include "dlopen_wrap.h"
#include "get_callable.h"
#include "hist_draw.h"
#include "wildcard_expand.h"
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <TChain.h>
#include <dlfcn.h>
#include <memory>
// #include <optional>
// #include <ranges>
#include <tuple>
// #include <unordered_map>

std::tuple<std::unique_ptr<TChain>, std::vector<std::unique_ptr<TChain>>>
prepare_chain(nlohmann::json &j) {
  std::tuple<std::unique_ptr<TChain>, std::vector<std::unique_ptr<TChain>>>
      ret{};
  std::vector<std::string> file_paths{};
  for (auto &file : j["files"]) {
    auto expanded_paths = expand_wildcard_path(file.get<std::string>());
    file_paths.insert(file_paths.end(), expanded_paths.begin(),
                      expanded_paths.end());
  }
  auto &&[filechain, friend_chains] = ret;
  filechain = std::make_unique<TChain>();
  for (auto &file : file_paths) {
    filechain->AddFile(
        (file + "?#" + j["treename"].get<std::string>()).c_str());
  }
  for (auto &friend_tree : j["friend_trees"]) {
    auto &chain = friend_chains.emplace_back(std::make_unique<TChain>());
    for (auto &file : file_paths) {
      chain->AddFile(
          (file + "?#" + friend_tree["treename"].get<std::string>()).c_str());
    }
    filechain->AddFriend(chain.get(), friend_tree.value("alias", "").c_str());
  }
  return ret;
}

std::vector<axis> get_bins(nlohmann::json &j) {
  std::vector<axis> ret{};
  if (j.find("bins_list") != j.end()) {
    for (auto &&bin_obj : j["bins_list"]) {
      ret.push_back(
          {bin_obj["var"], bin_obj["min"], bin_obj["max"], bin_obj["bins"]});
    }
    return ret;
  } else {
    if (j.find("var") != j.end() && j.find("bins") != j.end()) { // 1D case
      ret.push_back({j["var"], j["xmin"], j["xmax"], j["bins"]});
      return ret;
    } else { // 2-3D case
      ret.push_back({j["varx"], j["xmin"], j["xmax"], j["binsx"]});
      if (j.find("vary") != j.end()) {
        ret.push_back({j["vary"], j["ymin"], j["ymax"], j["binsy"]});
      }
      if (j.find("varz") != j.end()) {
        ret.push_back({j["varz"], j["zmin"], j["zmax"], j["binsz"]});
      }
      return ret;
    }
  }
}

auto analysis_entry_handle(ROOT::RDF::RNode preprocessed_node,
                           nlohmann::json &analysis_entry) {
  auto func = get_node_process_callable(analysis_entry["func"]);
  if (!func)
    throw std::runtime_error("Failed to get analysis function");
  auto result_node = (*func)(preprocessed_node);
  std::tuple<
      std::string, std::vector<ROOT::RDF::RResultPtr<TH1>>,
      std::vector<std::pair<THStack, std::vector<ROOT::RDF::RResultPtr<TH1D>>>>,
      // std::vector<ROOT::RDF::RResultPtr<TH2D>>,
      ROOT::RDF::RResultPtr<
          ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>>
      ret{};
  auto &&[plot_file, hist1ds, plots, snapshot_action] = ret;
  plot_file = analysis_entry["plot_file"].get<std::string>();

  auto cutset = analysis_entry.value("cutset", nlohmann::json{});
  std::vector<std::tuple<std::string, ROOT::RDF::RNode>> cut_nodes{};
  for (auto &&cuts : cutset) {
    std::vector<std::tuple<std::string, ROOT::RDF::RNode>>
        cut_nodes_next_layer{};
    for (auto &&cut : cuts) {
      std::cout << "cut: " << cut["cut"] << std::endl;

      cut_nodes_next_layer.emplace_back(cut["name"],
                                        do_cut(result_node, cut["cut"]));
      for (auto &&[name, node] : cut_nodes) {
        cut_nodes_next_layer.emplace_back(name + "_" +
                                              cut["name"].get<std::string>(),
                                          do_cut(node, cut["cut"]));
      }
    }
    cut_nodes.insert(cut_nodes.end(), cut_nodes_next_layer.begin(),
                     cut_nodes_next_layer.end());
  }

  // plot 1D hists
  for (auto &plotentry : analysis_entry["plots"]) {
    std::string name = plotentry["name"];
    auto bin_obj = get_bins(plotentry);
    auto cut = plotentry.value("cut", nlohmann::json{});
    std::string wname = plotentry.value("wname", "");
    hist1ds.emplace_back(draw_hists_nd(result_node, name, bin_obj, cut, wname));
    for (auto &[cutname, cutnode] : cut_nodes) {
      hist1ds.emplace_back(
          draw_hists_nd(cutnode, name + "_" + cutname, bin_obj, cut, wname));
    }
  }

  // Plot Hist2D
  for (auto &plotentry : analysis_entry["plot_2d"]) {
    std::string name = plotentry["name"];
    auto bin_obj = get_bins(plotentry);
    auto cut = plotentry.value("cut", nlohmann::json{});
    std::string wname = plotentry.value("wname", "");
    hist1ds.emplace_back(draw_hists_nd(result_node, name, bin_obj, cut, wname));
    for (auto &[cutname, cutnode] : cut_nodes) {
      hist1ds.emplace_back(
          draw_hists_nd(cutnode, name + "_" + cutname, bin_obj, cut, wname));
    }
  }

  // plot Stacks
  for (auto &plotentry : analysis_entry["stack_plots"]) {
    std::string var = plotentry["var"];
    std::string name = plotentry["name"];
    double xmin = plotentry.value("xmin", 0);
    double xmax = plotentry.value("xmax", 0);
    int nbins = plotentry.value("nbins", 128);
    std::string wname = plotentry.value("wname", "");
    auto &&[hs, histvec] = plots.emplace_back(
        std::make_pair(THStack{(name + "hs").c_str(), var.c_str()},
                       std::vector<ROOT::RDF::RResultPtr<TH1D>>{}));
    axis x{var, xmin, xmax, nbins};
    for (auto &cutentry : plotentry["cuts"]) {
      histvec.emplace_back(draw_hists(result_node, name, x, cutentry, wname));
    }
  }

  // for (auto &plot_2d_entry : analysis_entry["plot_3d"]) {
  //   std::string name = plot_2d_entry["name"];
  //   std::string varx = plot_2d_entry["varx"];
  //   std::string vary = plot_2d_entry["vary"];
  //   std::string varz = plot_2d_entry["varz"];
  //   double xmin = plot_2d_entry.value("xmin", 0.0);
  //   double xmax = plot_2d_entry.value("xmax", 0.0);
  //   double ymin = plot_2d_entry.value("ymin", 0.0);
  //   double ymax = plot_2d_entry.value("ymax", 0.0);
  //   double zmin = plot_2d_entry.value("zmin", 0.0);
  //   double zmax = plot_2d_entry.value("zmax", 0.0);
  //   int nbinsx = plot_2d_entry.value("nbinsx", 128);
  //   int nbinsy = plot_2d_entry.value("nbinsy", 128);
  //   int nbinsz = plot_2d_entry.value("nbinsz", 128);
  //   axis x{varx, xmin, xmax, nbinsx};
  //   axis y{vary, ymin, ymax, nbinsy};
  //   axis z{varz, zmin, zmax, nbinsz};
  //   auto cut = plot_2d_entry.value("cut", nlohmann::json{});
  //   std::string wname = plot_2d_entry.value("wname", "");
  //   hist1ds.emplace_back(draw_hists_3d(result_node, name, x, y, z, cut,
  //   wname)); for (auto &[cutname, cutnode] : cut_nodes) {
  //     hist1ds.emplace_back(
  //         draw_hists_3d(cutnode, name + "_" + cutname, x, y, z, cut, wname));
  //   }
  // }

  {
    auto &&collist = analysis_entry["output"].get<std::vector<std::string>>();
    if (!collist.empty()) {
      ROOT::RDF::RSnapshotOptions opts;
      opts.fLazy = true;
      opts.fOverwriteIfExists = true;
      snapshot_action = result_node.Snapshot(
          analysis_entry["treename"].get<std::string>(),
          analysis_entry["filename"].get<std::string>(), collist, opts);
    }
  }
  return ret;
}

void plugin_handle(ROOT::RDF::RNode rootnode, nlohmann::json &entry) {
  std::cout << "Processing " << entry["name"] << std::endl;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(rootnode);
#endif
  std::string so_name = entry["plugin"];
  dlopen_wrap so(so_name);

  std::unique_ptr<ProcessNodeI> preprocess =
      entry.find("preprocess") != entry.end()
          ? get_node_process_callable(entry["preprocess"])
          : std::make_unique<noop>();
  if (!preprocess)
    throw std::runtime_error("Failed to get preprocess function");

  auto preprocessed_node = (*preprocess)(rootnode);

  // double norm_factor{};
  std::unique_ptr<NormalizeI> normalize_func{};
  if (entry.find("normalize") != entry.end()) {
    auto normalize_conf = entry["normalize"];
    normalize_func = get_normalize_callable(normalize_conf);
    if (!normalize_func)
      throw std::runtime_error("Failed to get normalize function");
  }

  std::vector<std::tuple<
      std::string, std::vector<ROOT::RDF::RResultPtr<TH1>>,
      std::vector<std::pair<THStack, std::vector<ROOT::RDF::RResultPtr<TH1D>>>>,
      // std::vector<ROOT::RDF::RResultPtr<TH2D>>,
      ROOT::RDF::RResultPtr<
          ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>>>
      analysis_result_handles{};
  for (auto &analysis_entry : entry["analysis"]) {
    analysis_result_handles.emplace_back(
        analysis_entry_handle(preprocessed_node, analysis_entry));
  }

  // event loop should be triggered here if normalize_func is defined
  auto normalize_factor = normalize_func ? (*normalize_func)(rootnode) : 1.;
  if (normalize_factor != 1.) {
    for (auto &&[_, hist1ds, stacks, __] : analysis_result_handles) {
      for (auto &&hist : hist1ds) {
        hist->Scale(normalize_factor, "WIDTH");
      }
      for (auto &&[_, hists] : stacks) {
        for (auto &&hist : hists) {
          hist->Scale(normalize_factor, "WIDTH");
        }
      }
    }
  }

  for (auto &&[filename, hist1ds, stacks, snapshot_action] :
       analysis_result_handles) {
    std::cerr << "Writing to " << filename << std::endl;
    auto file = std::make_unique<TFile>(filename.c_str(), "RECREATE");
    file->cd();
    for (auto &&hist : hist1ds) {
      hist->Write();
    }
    for (auto &&[hs, hists] : stacks) {
      for (auto &&hist : hists) {
        hs.Add(hist.GetPtr());
      }
      hs.Write();
    }
    file->Write();
    file->Close();
  }
}
