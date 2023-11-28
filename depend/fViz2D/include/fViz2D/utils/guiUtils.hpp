#pragma once

#include "fViz2D/API/GUI_API.hpp"
#include "fViz2D/layers/imGuiLayer.hpp"

#include "fViz2D/renderer/rendererData.hpp"
#include "fViz2D/utils/imageUtils.hpp"

namespace GUI {

    void imGuiDrawTexture(TEX::textureID_t* texID_ptr, const char* windowName = "Texture Drawing");

    void imGuiCreateMenu(menuDefinition_t menuDefinition);

    menuDefinition_t getTestMenuDefinition(bool* testBool_ptr,  float* clearColorFirstElement_ptr, 
                                           float* noiseTintColorFirstElement_ptr);
}

