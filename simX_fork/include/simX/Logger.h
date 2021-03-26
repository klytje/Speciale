//
// Created by munk on 14-02-16.
//

#ifndef SIMX_LOGGER_H
#define SIMX_LOGGER_H

#include <spdlog/spdlog.h>

namespace simX {
    namespace log {
        using Level = spdlog::level::level_enum;
        using Logger = spdlog::logger;

        Level getLogLevel();
        void  setLogLevel(Level level);

        std::shared_ptr<Logger> getLogger(const std::string &name);
    }
}
#endif //SIMX_LOGGER_H
