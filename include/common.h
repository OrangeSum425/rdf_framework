#pragma once

#ifndef _COMMON_HEADER_
#define _COMMON_HEADER_
#include "ROOT/RDF/InterfaceUtils.hxx"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <functional>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

template <typename T> class ProcessNodeT {
public:
  virtual T operator()(ROOT::RDF::RNode) = 0;
  virtual void configure(const nlohmann::json &) {}
  ProcessNodeT() = default;
  ProcessNodeT(const ProcessNodeT &) = default;
  ProcessNodeT(ProcessNodeT &&) = default;
  ProcessNodeT &operator=(const ProcessNodeT &) = default;
  ProcessNodeT &operator=(ProcessNodeT &&) = default;
  virtual ~ProcessNodeT() = default;
};

template <typename T> class FactoryT {
public:
  std::unique_ptr<ProcessNodeT<T>> create(const std::string &name) const {
    auto it = creator.find(name);
    if (it == creator.end()) {
      return nullptr;
    }
    return it->second();
  }

  bool
  registerCreator(const std::string &name,
                  std::function<std::unique_ptr<ProcessNodeT<T>>()> m_creator) {
    this->creator[name] = m_creator;
    return true;
  }

  static FactoryT<T> &instance() {
    static FactoryT<T> instance;
    return instance;
  }

private:
  std::unordered_map<std::string,
                     std::function<std::unique_ptr<ProcessNodeT<T>>()>>
      creator;
  FactoryT() = default;
  FactoryT(const FactoryT &) = delete;
  FactoryT(FactoryT &&) = delete;
};

using ProcessNodeI = ProcessNodeT<ROOT::RDF::RNode>;
using ProcessNodeFactory = FactoryT<ROOT::RDF::RNode>;
using NormalizeI = ProcessNodeT<double>;
using NormalizeFactory = FactoryT<double>;

#define REGISTER_PROCESS_NODE(type)                                            \
  namespace {                                                                  \
  static bool type##registered =                                               \
      ProcessNodeFactory::instance().registerCreator(                          \
          #type, []() -> std::unique_ptr<ProcessNodeI> {                       \
            return std::make_unique<type>();                                   \
          });                                                                  \
  }

#define REGISTER_NORMALIZE(type)                                               \
  namespace {                                                                  \
  static bool type##registered = NormalizeFactory::instance().registerCreator( \
      #type, []() -> std::unique_ptr<NormalizeI> {                             \
        return std::make_unique<type>();                                       \
      });                                                                      \
  }

#endif