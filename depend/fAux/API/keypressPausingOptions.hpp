#pragma once

//This defines options for keypress pausing, including as a result of logging
//NOTE: changing these may lead to full recompile : (

//All LOG_CRITICALs and above will pause IF pauses are enabled (default 1)
#define PAUSE_ON_CRITICALS 1

//All LOG_ERRORs and above will pause IF pauses are enabled  (default 1)
#define PAUSE_ON_ERRORS 1

//The default for the following definitions is 0

//All LOG_WARNINGs and above will pause IF pauses are enabled
#define PAUSE_ON_WARNINGS 0

//All LOG_TRACE logs will be supressed if this is set to anything but 0
#define SUPRESS_TRACES 0

//Setting this to anything but 0 makes release supress only LOG_TRACEs
#define RELEASE_SUPRESS_ONLY_TRACES 0

//Setting this to anything but 0 makes release as verbose as debug
#define VERBOSE_RELEASE 0

//Setting this to anything but 0 makes debug NOT ASK for keypresses on GETCHAR_PAUSE
#define DONT_ASK_KEYPRESS_DEBUG 0

//Setting this to anything but 0 makes release ASK for keypresses on GETCHAR_PAUSE
#define ASK_KEYPRESS_ON_RELEASE 0

//Setting this to anything but 0 makes GETCHAR_FORCE_PAUSE act the same as GETCHAR_PAUSE
#define DONT_FORCE_KEYPRESS 0