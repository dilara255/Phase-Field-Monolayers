#pragma once

#include "fViz2D/API/GUI_API.hpp"
#include "fViz2D/API/returnCodes.hpp"

#include "fViz2D/API/imageUtils.hpp" //TODO: I don't want this to include any non-api header (I think?)

namespace F_V2 {

	//Returns BAD_DYNAMIC_DATA_FORMAT dynamicData_ptr's format is not supported, otherwise returns OK
	//Supported formats are the same as those supported by the texture loading functions in renderData.hpp
	F_V2::rendererRetCode_st rendererMain(IMG::generic2DfieldPtr_t* dynamicData_ptr,
                                          COLOR::rgbaF_t* clearColor_ptr,
		                                  GUI::menuDefinition_t userMenuDef = {},
		                                  GUI::filenameCallback_func* filenameFunc = nullptr,
                                          COLOR::colorInterpolation_t* scheme_ptr = nullptr,
                                          std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
                                          int width = 800, int height = 600,
                                          const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	//Sets returnCode_ptr to BAD_DYNAMIC_DATA_FORMAT if dynamicData_ptr's format is not supported, otherwise to OK
	//Supported formats are the same as those supported by the texture loading functions in renderData.hpp
	void rendererMainForSeparateThread(IMG::generic2DfieldPtr_t* dynamicData_ptr, 
		                               F_V2::rendererRetCode_st* returnCode_ptr, 
									   COLOR::rgbaF_t* clearColor_ptr,
		                               GUI::menuDefinition_t userMenuDef = {},
		                               GUI::filenameCallback_func* filenameFunc = nullptr,
                                       COLOR::colorInterpolation_t* scheme_ptr = nullptr,
                                       std::string windowName = "Ogl3 Render Test - imGui + Glfw", 
                                       int width = 800, int height = 600,
                                       const char* bannerPathFromBinary = F_V2::testBannerPathFromBinary);

	void saveCurrentImage();
}