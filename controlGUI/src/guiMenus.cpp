#include "imgui/imgui.h"

#include "fAux/API/miscStdHeaders.h"
#include "fViz2D/API/GUI_API.hpp"

#include "PFM_data.hpp"
#include "PFM_API.hpp"

#include "guiMenus.hpp"


static PFM::checkData_t* g_checkData_ptr = nullptr;
static const PFM::simConfig_t* g_simConfig_ptr = nullptr;
static PFM::simParameters_t* g_simParams_ptr = nullptr;

static PFM::simConfig_t* g_simConfigNextStart_ptr = nullptr;
static PFM::simParameters_t* g_simParamsNextStart_ptr = nullptr;

void checksMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_checkData_ptr == nullptr) { return; }

	ImGui::Text("%s", g_checkData_ptr->getChecksStr().c_str());
}

void controlFlowMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_simConfigNextStart_ptr == nullptr) { return; }
	if(g_simParamsNextStart_ptr == nullptr) { return; }
	if(g_simConfig_ptr == nullptr) { return; }
	if(g_simParams_ptr == nullptr) { return; }

	//TODO: make pause/resume and restart actually work (PFM/GUI API support needed)

	bool paused = PFM::isSimulationPaused();
	if (paused) {
		if (ImGui::Button("Resume")) { PFM::resumeSimulation(); }
	}
	else {
		if (ImGui::Button("Pause ")) { PFM::pauseSimulation(); }
	}
	ImGui::SameLine();
	if (ImGui::Button("Restart w/ these values")) { puts("\nrestart!\n"); }

	if (ImGui::Button("Copy Config")) {
		*g_simConfigNextStart_ptr = *g_simConfig_ptr;
		g_simConfigNextStart_ptr->stepsRan = 0;
	}
	ImGui::SameLine();
	if (ImGui::Button("Copy Params")) {
		*g_simParamsNextStart_ptr = *g_simParams_ptr;
	}

	//TODO: Get rid of the strings and instead expose the actual field for change

	ImGui::Text("Configuration for restart:\n%s", g_simConfigNextStart_ptr->getSimDataString().c_str());
	ImGui::Text("Parameters for restart:\n%s", g_simParamsNextStart_ptr->getSimParamsString().c_str());
}

void configAndParamsMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_simConfig_ptr == nullptr) { return; }
	if(g_simParams_ptr == nullptr) { return; }

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

GUI::menuDefinition_t PFM_GUI::getcontrolFlowMenuDefinition(PFM::simConfig_t* simConfigFromMain_ptr, 
		                                                    PFM::simParameters_t* simParamsFromMain_ptr) {
	g_simConfigNextStart_ptr = simConfigFromMain_ptr;
	g_simParamsNextStart_ptr = simParamsFromMain_ptr;
	return { controlFlowMenuFunc, "Control Flow menu" };
}