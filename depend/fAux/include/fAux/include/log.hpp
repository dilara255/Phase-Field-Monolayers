#pragma once

//Declaração funções log. Usa Spdlog. Exposto via macros em logAPI.hpp

#include "fAux/API/miscStdHeaders.h"

#include "fAux/API/core.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace az {

	class Log
	{
	public:
		static void init();
		inline static std::shared_ptr<spdlog::logger>& GetV2Logger() {
			if (!initialized) init();  return s_fViz2DLogger;
		}
		inline static std::shared_ptr<spdlog::logger>& GetCALogger() {
			if (!initialized) init(); return s_ClientAppLogger;
		}
		inline static std::shared_ptr<spdlog::logger>& GetDALogger() {
			if (!initialized) init(); return s_DebugAuxLogger;
		}
		inline static void log(std::shared_ptr<spdlog::logger> logger, const int degree,
			                   const char* file, const int line, 
							   const char* message, uint32_t trailingNewlines = 0);

	private:
		static std::shared_ptr<spdlog::logger> s_DebugAuxLogger; //used to provide extra info on Debug / Verbose
		static std::shared_ptr<spdlog::logger> s_fViz2DLogger;
		static std::shared_ptr<spdlog::logger> s_ClientAppLogger;		
		static int initialized;
	};

}

