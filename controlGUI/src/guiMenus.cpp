#include "imgui/imgui.h"

#include "fAux/API/miscStdHeaders.h"
#include "fViz2D/API/GUI_API.hpp"

#include "PFM_data.hpp"
#include "PFM_API.hpp"

#include "guiMenus.hpp"


static PFM::checkData_t* g_checkData_ptr = nullptr;
static const PFM::simConfig_t* g_simConfig_ptr = nullptr;
static PFM::simParameters_t* g_simParams_ptr = nullptr;

void checksMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_checkData_ptr == nullptr) { return; }

	ImGui::Text("%s", g_checkData_ptr->getChecksStr().c_str());
}

void configAndParamsMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {

	ImGui::Text("Simulation configuration:\n%s", g_simConfig_ptr->getSimDataString().c_str());
	ImGui::Text("Simulation parameters:\n%s", g_simParams_ptr->getSimParamsString().c_str());

	double completelyArbitraryMaxKplusATimesDt = 1.1;

	//These will be used to find the maximum allowed value of each parameter
	double minValuesForLimitCalc = 0.000001;
	double dividendForLimitCalc = std::max(minValuesForLimitCalc, g_simParams_ptr->dt);

	double mink = 0.0;

	double maxk = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc - g_simParams_ptr->A;
	ImGui::SliderScalar("k", ImGuiDataType_Double, &g_simParams_ptr->k, &mink, &maxk);

	double minA = 0.0;
	double maxA = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc - g_simParams_ptr->k;
	ImGui::SliderScalar("A", ImGuiDataType_Double, &g_simParams_ptr->A, &minA,&maxA);

	double minDt = 0.0;
	dividendForLimitCalc = std::max(minValuesForLimitCalc, g_simParams_ptr->A + g_simParams_ptr->k);
	double maxDt = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc;
	ImGui::SliderScalar("dt", ImGuiDataType_Double, &g_simParams_ptr->dt, &minDt,&maxDt);

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