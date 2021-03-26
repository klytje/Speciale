//
// Created by munk on 14-02-16.
//

#include "simX/Logger.h"

using namespace spdlog;
using namespace simX::log;

namespace {
    Level _level = spdlog::level::notice;
}

Level simX::log::getLogLevel() {
    return ::_level;
}

void ::simX::log::setLogLevel(Level lev) {
    _level = lev;
}

std::shared_ptr<Logger> simX::log::getLogger(const std::string &name) {
    auto old = spdlog::get(name);
    if (old) return old;
    
    auto log = spdlog::stderr_logger_st(name);
    log->set_level(_level);
    return log;
}
