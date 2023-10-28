#pragma once

/*
Define macros de API para facilitar log em console. Usa spdlog.
No momento não dá pra passar variáveis nem gravar em arquivo, 
apesar da biblioteca suportar isso.

TODO: Generalize the whole macro thingy and all that's associated with it
*/

#include "core.hpp"
#include "keypressPausingOptions.hpp"

//para os outros projetos poderem linkar as funções, declaradas em log.hpp
#ifndef F_AUX
    #include "../include/log.hpp"
#endif

namespace az {
    //Usage: LOG_*LEVEL*("message", int extraNewlinesAfter = 0) logs:
    //- Depending on level, file path and source line;
    //- Message, color coded by level and with time and project stamps; and
    //- Adds extraNewlinesAfter the message.
    void log(std::shared_ptr<spdlog::logger> logger, const int degree, const char* file, 
                                  const int line, const bool pause, const char* message, 
                                    uint32_t trailingNewlines = 0, bool supress = false);

    uint16_t RGB24toRGB565(uint8_t r, uint8_t g, uint8_t b);
}

namespace az {
    
}

#define L_TRACE 0
#define L_DEBUG 1
#define L_INFO 2
#define L_WARN 3
#define L_ERROR 4
#define L_CRITICAL 5

//Only macros from here on : )


//MACROS for user interaction

inline void consumeAllcharacters(bool showMessage = true, int trailingNewlines = 0) {
    int ch; //int to be able to tell apart EOF
    if (showMessage) {
        for(int i = 0; i < trailingNewlines; i++) {puts("\n"); }
        fputs("\t\tpress enter...", stdout);
    }
    while( (ch = getchar()) != '\n' && ch != EOF );
}

#if (defined AS_DEBUG && !DONT_ASK_KEYPRESS_DEBUG) || (defined AS_RELEASE && ASK_KEYPRESS_ON_RELEASE)
	#define GETCHAR_PAUSE consumeAllcharacters()
    #define GETCHAR_PAUSE_SILENT consumeAllcharacters(false)
#else
	#define GETCHAR_PAUSE puts("\n")
#endif // AS_DEBUG

#if (DONT_FORCE_KEYPRESS)
    #define GETCHAR_FORCE_PAUSE GETCHAR_PAUSE;
    #define GETCHAR_FORCE_PAUSE_SILENT GETCHAR_PAUSE_SILENT;
#else
	#define GETCHAR_FORCE_PAUSE consumeAllcharacters()
    #define GETCHAR_PAUSE_SILENT consumeAllcharacters(false)
#endif // AS_DEBUG

//MACROS log Agent System

#ifdef F_VIZ2D
    #define GETLOGGER az::Log::GetVZLogger()
#endif

#ifdef F_TESTAPP
    #define GETLOGGER az::Log::GetTALogger()
#endif

//Will only actually pause if GETCHAR_FORCE_PAUSE should pause
#define SHOULD_CRITICALS_PAUSE (PAUSE_ON_CRITICALS || PAUSE_ON_ERRORS || PAUSE_ON_WARNINGS)
#define SHOULD_ERRORS_PAUSE (PAUSE_ON_ERRORS || PAUSE_ON_WARNINGS)
#define SHOULD_WARNINGS_PAUSE PAUSE_ON_WARNINGS

//Usage: LOG_CRITICAL("message", int extraNewlinesAfter = 0)
#define LOG_CRITICAL(...) az::log(GETLOGGER, L_CRITICAL,\
                                   __FILE__, __LINE__, SHOULD_CRITICALS_PAUSE, __VA_ARGS__)

#ifdef AS_DISTRO
    #define LOG_ERROR(...)
#else
    //Usage: LOG_ERROR("message", int extraNewlinesAfter = 0)
    #define LOG_ERROR(...) az::log(GETLOGGER, L_ERROR,\
                                   __FILE__, __LINE__, SHOULD_ERRORS_PAUSE ,__VA_ARGS__)
#endif

#if !VERBOSE_RELEASE && defined AS_RELEASE
    #define LOG_WARN(...)
#else
    //Usage: LOG_WARN("message", int extraNewlinesAfter = 0)
    #define LOG_WARN(...) az::log(GETLOGGER, L_WARN,\
                                  __FILE__, __LINE__, SHOULD_WARNINGS_PAUSE, __VA_ARGS__)
#endif

#ifdef AS_DISTRO
    #define LOG_INFO(...)
#else
    //Usage: LOG_INFO("message", int extraNewlinesAfter = 0)
    #define LOG_INFO(...) az::log(GETLOGGER, L_INFO,\
                                  __FILE__, __LINE__, false, __VA_ARGS__)
#endif

#if (!VERBOSE_RELEASE && defined AS_RELEASE) || SUPRESS_TRACES
    #define LOG_TRACE(...)
#else
    //Usage: LOG_TRACE("message", int extraNewlinesAfter = 0)
    #define LOG_TRACE(...) az::log(GETLOGGER, L_TRACE,\
                                   __FILE__, __LINE__, false, __VA_ARGS__)
#endif

#if !(RELEASE_SUPRESS_ONLY_TRACES || VERBOSE_RELEASE) && defined AS_RELEASE
    #define LOG_DEBUG(...)
#else
    //Usage: LOG_TRACE("message", int extraNewlinesAfter = 0)
    #define LOG_DEBUG(...) az::log(GETLOGGER, L_DEBUG,\
                                   __FILE__, __LINE__, false, __VA_ARGS__)
#endif

