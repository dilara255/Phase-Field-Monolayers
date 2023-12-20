#pragma once

#include "fViz2D/prototypes.hpp"
#include "PFM_data.hpp"

namespace PFM_GUI {

	extern bool g_restartSimulationAfterStopped;
	extern bool g_dtLoweredForFirstSteps;
	extern double g_originalDt;

	GUI::menuDefinition_t getChecksMenuDefinition(PFM::checkData_t* checks_ptr);
	GUI::menuDefinition_t getConfigAndParamsMenuDefinition(const PFM::simConfig_t* simConfig_ptr, 
		                                                   PFM::simParameters_t* simParams_ptr);
	GUI::menuDefinition_t getcontrolFlowMenuDefinition(PFM::simConfig_t* simConfigFromMain_ptr, 
		                                               PFM::simParameters_t* simParamsFromMain_ptr);

	static const char* pfmGuiBannerPathFromBinary  = "../../../res/banners/bannerPFM.jpg";
}