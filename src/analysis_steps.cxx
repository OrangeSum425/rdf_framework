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
#include <tuple>
#include <unordered_map>

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

auto analysis_entry_handle(ROOT::RDF::RNode preprocessed_node,
                           nlohmann::json &analysis_entry) {
  auto func = get_node_process_callable(analysis_entry["func"]);
  if (!func)
    throw std::runtime_error("Failed to get analysis function");
  auto result_node = (*func)(preprocessed_node);
  std::tuple<
      std::string, std::vector<ROOT::RDF::RResultPtr<TH1D>>,
      std::vector<std::pair<THStack, std::vector<ROOT::RDF::RResultPtr<TH1D>>>>,
      std::vector<ROOT::RDF::RResultPtr<TH2D>>,
      ROOT::RDF::RResultPtr<
          ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>>
      ret{};
  auto &&[plot_file, hist1ds, plots, hist2ds, snapshot_action] = ret;
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
    std::string var = plotentry["var"];
    std::string name = plotentry["name"];
    double xmin = plotentry.value("xmin", 0.0);
    double xmax = plotentry.value("xmax", 0.0);
    int nbins = plotentry.value("nbins", 260);
    auto cut = plotentry.value("cut", nlohmann::json{});
    std::string wname = plotentry.value("wname", "");
    hist1ds.emplace_back(
        draw_hists(result_node, var, name, xmin, xmax, nbins, cut, wname));
    for (auto &[cutname, cutnode] : cut_nodes) {
      hist1ds.emplace_back(draw_hists(cutnode, var, name + "_" + cutname, xmin,
                                      xmax, nbins, cut, wname));
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
    for (auto &cutentry : plotentry["cuts"]) {
      histvec.emplace_back(draw_hists(result_node, var, name, xmin, xmax, nbins,
                                      cutentry, wname));
    }
  }

  // Plot Hist2D
  for (auto &plot_2d_entry : analysis_entry["plot_2d"]) {
    std::string varx = plot_2d_entry["varx"];
    std::string vary = plot_2d_entry["vary"];
    std::string name = plot_2d_entry["name"];
    double xmin = plot_2d_entry.value("xmin", 0.0);
    double xmax = plot_2d_entry.value("xmax", 0.0);
    double ymin = plot_2d_entry.value("ymin", 0.0);
    double ymax = plot_2d_entry.value("ymax", 0.0);
    int nbinsx = plot_2d_entry.value("nbinsx", 128);
    int nbinsy = plot_2d_entry.value("nbinsy", 128);
    auto cut = plot_2d_entry.value("cut", nlohmann::json{});
    std::string wname = plot_2d_entry.value("wname", "");
    hist2ds.emplace_back(draw_hists_2d(result_node, varx, vary, name, xmin,
                                       xmax, ymin, ymax, nbinsx, nbinsy, cut,
                                       wname));
    for (auto &[cutname, cutnode] : cut_nodes) {
      hist2ds.emplace_back(draw_hists_2d(cutnode, varx, vary,
                                         name + "_" + cutname, xmin, xmax, ymin,
                                         ymax, nbinsx, nbinsy, cut, wname));
    }
  }
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
      std::string, std::vector<ROOT::RDF::RResultPtr<TH1D>>,
      std::vector<std::pair<THStack, std::vector<ROOT::RDF::RResultPtr<TH1D>>>>,
      std::vector<ROOT::RDF::RResultPtr<TH2D>>,
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
    for (auto &&[_, hist1ds, stacks, hist2ds, __] : analysis_result_handles) {
      for (auto &&hist : hist1ds) {
        hist->Scale(normalize_factor, "WIDTH");
      }
      for (auto &&[_, hists] : stacks) {
        for (auto &&hist : hists) {
          hist->Scale(normalize_factor, "WIDTH");
        }
      }
      for (auto &&hist : hist2ds) {
        hist->Scale(normalize_factor, "WIDTH");
      }
    }
  }

  for (auto &&[filename, hist1ds, stacks, hist2ds, snapshot_action] :
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
    for (auto &&hist : hist2ds) {
      hist->Write();
    }
    file->Write();
    file->Close();
  }
}
