#include "dlopen_wrap.h"
#include <iostream>

dlopen_wrap::dlopen_wrap(const std::string& filename, int flag) : handle(dlopen(filename.c_str(), flag)) {
  if (!handle) {
    std::cerr << "Cannot open library: " << dlerror() << '\n';
    throw std::runtime_error("Cannot open library");
  }
}

dlopen_wrap::~dlopen_wrap() {
  if (handle)
    dlclose(handle);
}

void *dlopen_wrap::get_handle() { return handle; }

void *dlopen_wrap::get_symbol(const char *symbol) {
  if (!handle)
    return nullptr;
  return dlsym(handle, symbol);
}

dlopen_wrap::dlopen_wrap(dlopen_wrap &&other)  noexcept : handle(other.handle) {
  other.handle = nullptr;
}

dlopen_wrap &dlopen_wrap::operator=(dlopen_wrap &&other) noexcept {
  handle = other.handle;
  other.handle = nullptr;
  return *this;
}