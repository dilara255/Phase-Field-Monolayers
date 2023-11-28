#pragma once

#include "fAux/API/core.hpp"

#include "fViz2D/utils/imageUtils.hpp"
#include "fViz2D/utils/guiUtils.hpp"
#include "fViz2D/renderer/rendererData.hpp"

#include "fViz2D/resourcePaths.hpp"

namespace F_V2 {

	//dynamicData_ptr should hold either a rgbaImage_t or a floats2Dfield_t (see imageUtils.hpp)
	//In case of a bad dynamicData_ptr type, returns BAD_DYNAMIC_DATA_FORMAT, otherwise returns OK
	//TODO: ADD TEST
	F_V2_API F_V2::rendererRetCode_st spawnRendererOnThisThread(IMG::generic2DfieldPtr_t* dynamicData_ptr,
												  COLOR::rgbaF_t* clearColor_ptr,
												  GUI::menuDefinition_t userMenuDef = {},
		                                          GUI::filenameCallback_func* filenameFunc = nullptr,
												  COLOR::colorInterpolation_t* scheme_ptr = nullptr,
												  std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
												  int width = 800, int height = 600,
												  const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	//dynamicData_ptr should hold either a rgbaImage_t or a floats2Dfield_t (see imageUtils.hpp)
	//In case of a bad dynamicData_ptr type, returns "empty" thread and sets returnCode to BAD_DYNAMIC_DATA_FORMAT
	F_V2_API std::thread spawnRendererOnNewThread(IMG::generic2DfieldPtr_t* dynamicData_ptr, 
											      F_V2::rendererRetCode_st* returnCode_ptr, 
											      COLOR::rgbaF_t* clearColor_ptr,
											      GUI::menuDefinition_t userMenuDef = {},
		                                          GUI::filenameCallback_func* filenameFunc = nullptr,
											      COLOR::colorInterpolation_t* scheme_ptr = nullptr,
											      std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
											      int width = 800, int height = 600,
											      const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	//TODO: add function to tell renderer to save current image
}