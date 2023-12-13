#pragma once

#include "fAux/API/core.hpp"
#include "fAux/API/miscStdHeaders.h"

#include "fViz2D/prototypes.hpp"

#include "fViz2D/utils/imageUtils.hpp"
#include "fViz2D/renderer/rendererData.hpp"

#include "fViz2D/resourcePaths.hpp"

namespace F_V2 {

	typedef void synchCallback_func();

	//dynamicData_ptr should hold either a rgbaImage_t or a floats2Dfield_t (see imageUtils.hpp)
	//In case of a bad dynamicData_ptr type, returns BAD_DYNAMIC_DATA_FORMAT, otherwise returns OK
	//TODO: ADD TEST
	F_V2_API F_V2::rendererRetCode_st spawnRendererOnThisThread(IMG::generic2DfieldPtr_t* dynamicData_ptr,
												  COLOR::rgbaF_t* clearColor_ptr,
												  GUI::menuDefinitionList_t* userMenuDefs_ptr = nullptr,
		                                          GUI::filenameCallback_func* filenameFunc = nullptr,
												  COLOR::colorInterpolation_t* scheme_ptr = nullptr,
												  std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
												  int width = 800, int height = 600,
											      bool createDefaultRendererMenu = true,
		                                          synchCallback_func synchCallback = defaultSynchCallback,
												  const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	//dynamicData_ptr should hold either a rgbaImage_t or a floats2Dfield_t (see imageUtils.hpp)
	//In case of a bad dynamicData_ptr type, returns "empty" thread and sets returnCode to BAD_DYNAMIC_DATA_FORMAT
	F_V2_API std::thread spawnRendererOnNewThread(IMG::generic2DfieldPtr_t* dynamicData_ptr, 
											      F_V2::rendererRetCode_st* returnCode_ptr, 
											      COLOR::rgbaF_t* clearColor_ptr,
											      GUI::menuDefinitionList_t* userMenuDefs_ptr = nullptr,
		                                          GUI::filenameCallback_func* filenameFunc = nullptr,
											      COLOR::colorInterpolation_t* scheme_ptr = nullptr,
											      std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
											      int width = 800, int height = 600,
		                                          bool createDefaultRendererMenu = true,
		                                          synchCallback_func synchCallback = defaultSynchCallback,
											      const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	//Used to expose some control over the renderer (eg, to create GUIs)
	typedef struct rendererControlPtrs_st {
		bool initialized = false;
		bool* shouldInterpolateColors_ptr = nullptr; 
		bool* keepRendering_ptr = nullptr; 
		void* ImGuiIO_ptr = nullptr; 
		bool* shouldSave_ptr = nullptr;
		bool* saveCalledFromGUI_ptr = nullptr;
		int* steps_ptr = nullptr;

		void loadPointers(bool* p_shouldInterpolateColors_ptr, bool* p_keepRendering_ptr, void* p_ImGuiIO_ptr, 
			              bool* p_saveCalledFromGUI_ptr, bool* p_shouldSave_ptr, int* p_steps_ptr) {

			shouldInterpolateColors_ptr = p_shouldInterpolateColors_ptr;
			keepRendering_ptr = p_keepRendering_ptr;
			ImGuiIO_ptr = p_ImGuiIO_ptr;
			shouldSave_ptr = p_shouldSave_ptr;
			saveCalledFromGUI_ptr = p_saveCalledFromGUI_ptr;
			steps_ptr = p_steps_ptr;
			initialized = true;
		}
	} rendererControlPtrs_t;

	//Quality will be clamped tp [0, 100]
	//Filename will use the naming function passsed to the renderer (or default), and tell it this is an API call
	void saveCurrentImage(int quality = 100);
}