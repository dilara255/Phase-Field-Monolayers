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

void createParameterSliders(int stepsRan, PFM::simParameters_t* param_ptr, double* GUIissuedDt_ptr) {
	
	ImGui::Checkbox("Adaptative dt", &param_ptr->useAdaptativeDt);

	PFM::parameterBounds_t bounds = 
		PFM::calculateParameterBounds(param_ptr->k, param_ptr->A, param_ptr->dt, stepsRan);

	if (param_ptr->useAdaptativeDt) {
		bounds.maxA = 1.0;
		bounds.maxK = 1.0;
	}

	ImGui::SliderScalar("k", ImGuiDataType_Double, &param_ptr->k, &bounds.minK, &bounds.maxK);
	ImGui::SliderScalar("A", ImGuiDataType_Double, &param_ptr->A, &bounds.minA, &bounds.maxA);
	if (!param_ptr->useAdaptativeDt) { 
		ImGui::SliderScalar("dt", ImGuiDataType_Double, GUIissuedDt_ptr, &bounds.minDt, &bounds.maxDt); 
	}
}

void controlFlowMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_simConfigNextStart_ptr == nullptr) { return; }
	if(g_simParamsNextStart_ptr == nullptr) { return; }
	if(g_simConfig_ptr == nullptr) { return; }
	if(g_simParams_ptr == nullptr) { return; }

	bool paused = PFM::isSimulationPaused();
	if (paused) {
		if (ImGui::Button("Resume")) { PFM::resumeSimulation(); }
	}
	else {
		if (ImGui::Button("Pause ")) { PFM::pauseSimulation(); }
	}
	ImGui::SameLine();
	ImGui::Checkbox("Restart after exit", &PFM_GUI::g_restartSimulationAfterStopped);

	if (ImGui::Button("Copy Config")) {
		*g_simConfigNextStart_ptr = *g_simConfig_ptr;
		g_simConfigNextStart_ptr->stepsRan = 0;
	}
	ImGui::SameLine();
	if (ImGui::Button("Copy Params")) {
		*g_simParamsNextStart_ptr = *g_simParams_ptr;
	}

	ImGui::Text("Configuration for restart:\n");
	//TODO-CRITICAL: Extract *A LOT* of this : )
	int minDimension = 1;
	int maxDimension = 1024;
	int minCells = 0;
	int maxCells = (maxDimension / 16) * (maxDimension / 16);
	double minBias = -0.5;
	double maxBias = 0.5;
	int minInitialCond = 0;
	int maxInitialCond = (int)PFM::initialConditions::TOTAL_INITIAL_CONDS - 1;
	int minMethod = 0;
	int maxMethod = (int)PFM::integrationMethods::TOTAL_METHODS - 1;
	//TODO: Change when the others are (re-)implemented
	int minSimFunc = (int)PFM::simFuncEnum::SINGLE_LAYER_CH_SIM;
	int maxSimFunc = (int)PFM::simFuncEnum::SINGLE_LAYER_CH_SIM;
	
	ImGui::InputInt("Width", &g_simConfigNextStart_ptr->width, 1, 64);
	g_simConfigNextStart_ptr->width = std::clamp(g_simConfigNextStart_ptr->width, minDimension, maxDimension);
	
	ImGui::InputInt("Height", &g_simConfigNextStart_ptr->height, 1, 64);
	g_simConfigNextStart_ptr->height = std::clamp(g_simConfigNextStart_ptr->height, minDimension, maxDimension);
	
	ImGui::InputInt("Cells", &g_simConfigNextStart_ptr->cells, 1, 8);
	g_simConfigNextStart_ptr->cells = std::clamp(g_simConfigNextStart_ptr->cells, minCells, maxCells);
	
	ImGui::SliderScalar("Bias", ImGuiDataType_Double, &g_simConfigNextStart_ptr->bias, &minBias, &maxBias);
	ImGui::SliderScalar("Initial Condition", ImGuiDataType_S32, &g_simConfigNextStart_ptr->initialContidion, 
		                                                                   &minInitialCond, &maxInitialCond);
	if (g_simConfigNextStart_ptr->initialContidion != PFM::initialConditions::LINEAR_RANDOM) {
		g_simConfigNextStart_ptr->cells = std::max(1, g_simConfigNextStart_ptr->cells);
	}
	ImGui::SliderScalar("Method", ImGuiDataType_S32, &g_simConfigNextStart_ptr->method, &minMethod, &maxMethod);
	ImGui::SliderScalar("Simulation", ImGuiDataType_S32, &g_simConfigNextStart_ptr->simulFunc, 
		                                                             &minSimFunc, &maxSimFunc);
	g_simConfigNextStart_ptr->perCellLayer = false;
	if(g_simConfigNextStart_ptr->simulFunc == PFM::simFuncEnum::MULTI_LAYER_CH_SIM) {
		g_simConfigNextStart_ptr->perCellLayer = true;
	}
	
	ImGui::Checkbox("Start Paused", &g_simConfigNextStart_ptr->startPaused);
	
	ImGui::Text("Parameters for restart:\n");
	
	createParameterSliders(2 * PFM::completelyArbitraryStepToUnlockFullDt, 
		                   g_simParamsNextStart_ptr, &(g_simParamsNextStart_ptr->dt));
	g_simConfigNextStart_ptr->stepsRan = 0; //change place

	g_simParamsNextStart_ptr->lambda = 
		PFM::getLambdaFromKandA(g_simParamsNextStart_ptr->k, g_simParamsNextStart_ptr->A);
	g_simParamsNextStart_ptr->gamma = 
		PFM::getGammaFromKandA(g_simParamsNextStart_ptr->k, g_simParamsNextStart_ptr->A);

	ImGui::Text("Lambda: %f\nGamma: %f", g_simParamsNextStart_ptr->lambda, g_simParamsNextStart_ptr->gamma);
}

void configAndParamsMenuFunc(F_V2::rendererControlPtrs_t* rendererCtrl_ptrs) {
	if(g_simConfig_ptr == nullptr) { return; }
	if(g_simParams_ptr == nullptr) { return; }

	ImGui::Text("Simulation configuration:\n%s", g_simConfig_ptr->getSimDataString().c_str());
	ImGui::Text("Simulation parameters:\n%s", g_simParams_ptr->getSimParamsString().c_str());

	static double GUIissuedDt = PFM::getIntendedDt();

	createParameterSliders(g_simConfig_ptr->stepsRan, g_simParams_ptr, &GUIissuedDt);

	if (GUIissuedDt != PFM::getIntendedDt()) {
		//A change was actually issued by the user, so:
		PFM::setIntendedDt(GUIissuedDt);
	}

	PFM::updatePhysicalParameters();

	double changeToSave = PFM::getMaxTotalChangePerElementPerCheckAdded();
	ImGui::InputDouble("Total Change/Element to save", &changeToSave, 0.0001, 0.005);
	changeToSave = std::clamp(changeToSave, 0.0000000000001, 0.5);
	PFM::setMaxTotalChangePerElementPerCheckAdded(changeToSave);

	int stepsPerCheck = PFM::getMaxStepsPerCheckAdded();
	ImGui::InputInt("Max steps/check", &stepsPerCheck, 10, 100);
	stepsPerCheck = std::clamp(stepsPerCheck, 1, INT32_MAX);
	PFM::setMaxStepsPerCheckAdded((uint32_t)stepsPerCheck);
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