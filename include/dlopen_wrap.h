#pragma once

#include <dlfcn.h>
#include <string_view>

class dlopen_wrap {
public:
  dlopen_wrap(std::string_view name, int flag = RTLD_NOW);
  ~dlopen_wrap();
  void *get_handle();
  void *get_symbol(const char *symbol);
  dlopen_wrap(dlopen_wrap &&) = default;
  dlopen_wrap &operator=(dlopen_wrap &&) = default;

private:
  void *handle{};
  dlopen_wrap(const dlopen_wrap &other) = delete;
  dlopen_wrap &operator=(const dlopen_wrap &other) = delete;
};