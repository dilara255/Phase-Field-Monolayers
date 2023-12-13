#pragma once

#include "fAux/API/miscStdHeaders.h"

namespace GUI { 
	typedef struct menuDefinition_st menuDefinition_t; 
	typedef std::string filenameCallback_func(int steps, bool calledFromGui);
	typedef std::vector<GUI::menuDefinition_t> menuDefinitionList_t;
}

namespace F_V2 {
	void defaultSynchCallback();
}