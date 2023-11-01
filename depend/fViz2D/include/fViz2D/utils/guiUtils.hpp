#pragma once

#include "fViz2D/layers/imGuiLayer.hpp"

#include "fViz2D/renderer/rendererData.hpp"
#include "fViz2D/utils/imageUtils.hpp"

namespace GUI {
	void imGuiDrawTexture(TEX::textureID_t* texID_ptr, const char* windowName = "Texture Drawing");

	void imGuiTestMenu(ImGuiIO& io, float* clearColorFirstElement_ptr, float* noiseTintColorFirstElement_ptr,
	                                                             bool* keepRendering_ptr, bool* testBool_ptr);
}

