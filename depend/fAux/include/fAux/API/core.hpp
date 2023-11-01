#pragma once

//Definições para controle dos projetos, sistemas, compilação, linkagem, etc
//TODO: Move _API definitions to the build system

#ifndef F_AUX
	#include "logAPI.hpp" 
#endif

#define BIT(x) (1 << x) //TODO: find a better place for this

namespace AZ{
	static long long getExpectedWakeUpDelay(long long microsToSleep = 0){
		//return (long long)round(478 - 0.00373*microsToSleep);
		return 0;
	}

	static const char* resPathFromBinary = "../../../res/";
}

///TODO: MOVE TO PLATFORM-LAYER (and maybe include that here):
#ifdef F_OS_WINDOWS

	#define AZ_LOG_TRACE_COLOR (FOREGROUND_BLUE | FOREGROUND_GREEN)
	#define AZ_LOG_DEBUG_COLOR (FOREGROUND_BLUE | FOREGROUND_RED | FOREGROUND_INTENSITY)
	#define AZ_LOG_INFO_COLOR (FOREGROUND_GREEN | FOREGROUND_INTENSITY)

#elif defined(F_OS_LINUX)

	const char azLogTraceColor[] = "\033[36m";
	const char azLogDebugColor[] = "\033[35m\033[1m";
	const char azLogInfoColor[] = "\033[32m\033[1m";

	#define AZ_LOG_TRACE_COLOR (azLogTraceColor)
	#define AZ_LOG_DEBUG_COLOR (azLogDebugColor)
	#define AZ_LOG_INFO_COLOR (azLogInfoColor)

#else
	#error "OS not recognized"
#endif
