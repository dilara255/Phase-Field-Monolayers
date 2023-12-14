#include "guiMenus.hpp"
#include "fViz2D/API/GUI_API.hpp"
#include "PFM_data.hpp"
#include "PFM_API.hpp"
#include "imgui/imgui.h"

static PFM::checkData_t* g_checkData_ptr = nullptr;
static const PFM::simConfig_t* g_simConfig_ptr = nullptr;
static PFM::simParameters_t* g_simParams_ptr = nullptr;

void checksMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	ImGui::Text("Checks @ step %d:\nDensity: %f\nchange: %f, absolute change: %f",
		        g_checkData_ptr->step, g_checkData_ptr->lastDensity, 
		        g_checkData_ptr->lastDensityChange, g_checkData_ptr->lastAbsoluteChange);
}

void configAndParamsMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {

	ImGui::Text("Simulation configuration:\n%s", g_simConfig_ptr->getSimDataString().c_str());
	ImGui::Text("Simulation parameters:\n%s", g_simParams_ptr->getSimParamsString().c_str());

	double completelyArbitraryMaxKplusA = 1.1;

	double mink = 0.0;
	double maxk = completelyArbitraryMaxKplusA - g_simParams_ptr->A;
	ImGui::SliderScalar("k", ImGuiDataType_Double, &g_simParams_ptr->k, &mink, &maxk);

	double minA = 0.0;
	double maxA = completelyArbitraryMaxKplusA - g_simParams_ptr->k;
	ImGui::SliderScalar("A", ImGuiDataType_Double, &g_simParams_ptr->A, &minA,&maxA);

	PFM::updatePhysicalParameters();
}

GUI::menuDefinition_t PFM_GUI::getChecksMenuDefinition(PFM::checkData_t* checks_ptr) {
	g_checkData_ptr = checks_ptr;
	return { checksMenuFunc, "Checks menu" };
}

GUI::menuDefinition_t PFM_GUI::getConfigAndParamsMenuDefinition(const PFM::simConfig_t* simConfig_ptr, 
		                                                        PFM::simParameters_t* simParams_ptr) {
	g_simConfig_ptr = simConfig_ptr;
	g_simParams_ptr = simParams_ptr;
	return { configAndParamsMenuFunc, "Config menu" };
}