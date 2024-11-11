#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <RVersion.h>
#include <TChain.h>
#include <common.h>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

#include "analysis_steps.h"
#include "env_handler.h"

int main(int argc, char **argv) {
  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <json_file>" << std::endl;
    return 1;
  }

  if (argc > 2) {
    for (int i = 2; i < argc; i++) {
      std::string arg{argv[i]};
      auto pos = arg.find('=');
      if (pos != std::string_view::npos) {
        env_handler::get_instance().set_env(arg.substr(0, pos),
                                            arg.substr(pos + 1));
      }
    }
  }

  nlohmann::json j;
  {
    std::ifstream json_file(argv[1]);
    std::stringstream buffer;
    buffer << json_file.rdbuf();
    std::stringstream{buffer.str() | env_handler::get_instance()} >> j;
  }
  // friend_chains is not actually used
  // but we need to keep it alive
  auto &&[filechain, friend_chains] = prepare_chain(j);
  ROOT::RDataFrame df(*filechain);

  auto rootnode = ROOT::RDF::AsRNode(df);

  auto &&plugins_entry = j["plugins"];
  if (plugins_entry.is_array()) {
      for (auto &&p : plugins_entry) {
          try {
              std::cout << "Attempting to run plugin: " << p["name"] << std::endl;
              plugin_handle(rootnode, p);
          } catch (std::exception &e) {
              std::cerr << "Failed to run plugin: " << p["name"].get<std::string>()
                        << "\n\t" << e.what() << "\n\t"
                        << "Skipping..." << std::endl;
          } catch (...) {
             std::cerr << "Failed to run plugin: " << p["name"].get<std::string>()
                       << "\n\t"
                       << "Skipping..." << std::endl;
          }
      }
  } else {
      try {
          std::cout << "Attempting to run plugin: " << plugins_entry["name"] << std::endl;
          plugin_handle(rootnode, plugins_entry);
      } catch (std::exception &e) {
         std::cerr << "Failed to run plugin: "
                    << plugins_entry["name"].get<std::string>() << "\n\t"
                    << e.what() << "\n\t"
                    << "Skipping..." << std::endl;
      } catch (...) {
          std::cerr << "Failed to run plugin: "
                    << plugins_entry["name"].get<std::string>() << "\n\t"
                    << "Skipping..." << std::endl;
      }
  }
  return 0;
}
