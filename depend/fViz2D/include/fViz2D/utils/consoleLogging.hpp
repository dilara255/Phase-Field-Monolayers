#pragma once

// USE fAux Logger:
#include "fAux/API/logAPI.hpp"

/* DON'T USE fAux Logger:

//Ignore this for now : )
//TODO: don't

//TODO: option parity with AZ_log and extract this to Aux (as header-only)
    
#include <cstdio>

#define LOG_CONSOLE(...) printf("%s\n", __VA_ARGS__)

#define LOG_CRITICAL(...) LOG_CONSOLE(...)
#define LOG_ERROR(...) LOG_CONSOLE(...)

#if defined AS_RELEASE
    #define LOG_WARN(...)
    #define LOG_DEBUG(...)
    #define LOG_TRACE(...)
#else
    #define LOG_WARN(...) LOG_CONSOLE(...)
    #define LOG_DEBUG(...) LOG_CONSOLE(...)
    #define LOG_TRACE(...) LOG_CONSOLE(...)
#endif

#define LOG_INFO(...) LOG_CONSOLE(...)

*/