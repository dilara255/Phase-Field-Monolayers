#pragma once

#include "fAux/API/miscStdHeaders.h"

namespace GUI {

	typedef std::vector<void*> hookList_t;
    typedef void guiMenu_func(hookList_t hooks);

    //TODO: add menu name? Note that imGui expects names to not be empty
    typedef struct menuDefinition_st {
        hookList_t hooks;
        guiMenu_func* menuFunc_ptr = nullptr;
        std::string menuName = "Client Menu";
    } menuDefinition_t;

    typedef std::string filenameCallback_func(int steps);
}