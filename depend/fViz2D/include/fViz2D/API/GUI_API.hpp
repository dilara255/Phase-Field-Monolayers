#pragma once

#include "fViz2D/API/FV2_API.hpp"

namespace GUI {

	typedef std::vector<void*> hookList_t;

	typedef void guiMenu_func(F_V2::rendererControlPtrs_t* rendererControl_ptr);

    typedef struct menuDefinition_st {
        guiMenu_func* menuFunc_ptr = nullptr;
        std::string menuName = "Client Menu"; //Note: Imgui expects name to not be empty. Are we handling that?
    } menuDefinition_t;

	menuDefinition_t getTestMenuDefinition(bool* testBool_ptr,  float* clearColorFirstElement_ptr, 
                                           float* noiseTintColorFirstElement_ptr);
}