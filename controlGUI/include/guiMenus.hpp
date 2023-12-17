#pragma once

#include "fViz2D/prototypes.hpp"
#include "PFM_data.hpp"

namespace PFM_GUI {
	GUI::menuDefinition_t getChecksMenuDefinition(PFM::checkData_t* checks_ptr);
	GUI::menuDefinition_t getConfigAndParamsMenuDefinition(const PFM::simConfig_t* simConfig_ptr, 
		                                                   PFM::simParameters_t* simParams_ptr);

	static const char* pfmGuiBannerPathFromBinary  = "../../../res/banners/bannerPFM.jpg";
}