#include "fViz2D/API/FV2_API.hpp"
#include "fViz2D/API/GUI_API.hpp"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"

#include "guiTests.hpp"
#include "guiControlMain.hpp"

#include "guiMenus.hpp"

#define RUN_GUI_TESTS 0
#define RUN_SIM_DATACTRL_TESTS 0
#define RUN_SINGLE_LAYER_CH_SIM 1
#define RUN_MULTI_LAYER_CH_SIM 0
#define TOTAL_SIM_FUNCS ((int)PFM::simFuncEnum::TOTAL_SIM_FUNCS)

//TODO: These should exist in the simulation project, not here
typedef struct simConfig_st {
	int width, height, cells;
	PFM::initialConditions initialCond = PFM::initialConditions::EVENLY_SPACED_INDEX;
	double bias = 0;
	PFM::integrationMethods method = PFM::integrationMethods::FTCS;
	double dt = 1;
	double lambda = 3; //3 //7.824813
	double gamma = 0.06; //0.06 //0.043986
} simConfig_t;

//TODO: These should exist in the simulation project, not here
simConfig_t defaultParamsPerSimulType[TOTAL_SIM_FUNCS] = {
	{512, 512, 50},
	{128, 128, 50, PFM::initialConditions::LINEAR_RANDOM, -0.3, PFM::integrationMethods::FTCS},
	{128, 128, 5}
};

int main() {

	bool result = true;

	if(RUN_GUI_TESTS) { result &= PFM_GUI_TESTS::guiLinkingAndDependencyTests(); }
	if(RUN_SIM_DATACTRL_TESTS) { result &= PFM_GUI::runSimulationWithGUI(PFM::simFuncEnum::DATA_CONTROL_TEST); }
	if(RUN_SINGLE_LAYER_CH_SIM) { result &= PFM_GUI::runSimulationWithGUI(PFM::simFuncEnum::SINGLE_LAYER_CH_SIM); }
	if(RUN_MULTI_LAYER_CH_SIM) { result &= PFM_GUI::runSimulationWithGUI(PFM::simFuncEnum::MULTI_LAYER_CH_SIM); }

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

//TODO: receive config and param structs
//Should also receive the number of steps to run and/or minimumAbsoluteChange
//Should receive stepsPerSave and/or totalAbsoluteChangePerSave as well
//Make sure if a posite number of steps was passed, manually closing the renderer still stops the simulation
//Add option for data saves to also save the image from the GUI
bool PFM_GUI::runSimulationWithGUI(PFM::simFuncEnum simulationFunctionToRun) {
	if ((int)simulationFunctionToRun >= TOTAL_SIM_FUNCS) {
		LOG_ERROR("Invalid simulation function selected");
		return false;
	}

	LOG_DEBUG("Will run the simulation and display on the GUI");

	//***********vvvv This is setting up and running the simulation vvvv*******************
	if(simulationFunctionToRun != PFM::simFuncEnum::DATA_CONTROL_TEST) { 
		LOG_WARN("The menu to validate this test is not supported on the current version of the controller. Will say that the manual inspection failed. TODO: Reimplement support");
	}
	int simIndex = (int)simulationFunctionToRun;
	
	//These should be received (or set to their defaults otherwise)
	int width = defaultParamsPerSimulType[simIndex].width;
	int height = defaultParamsPerSimulType[simIndex].height;
	int cells = defaultParamsPerSimulType[simIndex].cells;
	auto initialCondition = defaultParamsPerSimulType[simIndex].initialCond;
	double bias = defaultParamsPerSimulType[simIndex].bias;	
	PFM::integrationMethods method = defaultParamsPerSimulType[simIndex].method;
	double dt = defaultParamsPerSimulType[simIndex].dt;
	double lambda = defaultParamsPerSimulType[simIndex].lambda;
	double gamma = defaultParamsPerSimulType[simIndex].gamma;

	PFM::fieldDimensions_t dimensions = {(size_t)width, (size_t)height};
	
	PFM::initializeSimulation(dimensions, cells, initialCondition, bias);
	LOG_INFO("Simulation initialized. Will run the simulation");
	PFM::runForSteps(-1, lambda, gamma, dt, simulationFunctionToRun, method);

	//***********^^^^ This is setting up and running the simulation ^^^^*******************

	//This should be received
	GUI::filenameCallback_func* filenameFunc = PFM::getFileName;

	//Forming the menus:
	//The menuList should be received 
	auto checks_ptr = PFM::getCheckDataPtr();
	auto config_ptr = PFM::getSimConfigPtr();
	auto params_ptr = PFM::getSimParamsPtr();

	if(checks_ptr == nullptr || config_ptr == nullptr || params_ptr == nullptr) {
		LOG_ERROR("Could get the pointers to form the menus");
		return false;
	}

	GUI::menuDefinitionList_t menuList;
	menuList.push_back(PFM_GUI::getChecksMenuDefinition(checks_ptr));
	menuList.push_back(PFM_GUI::getConfigAndParamsMenuDefinition(config_ptr, params_ptr));

	//Loading the color scheme:
	//This should also be received
	COLOR::colorInterpolation_t scheme;
	scheme.loadScheme(&COLOR::defaultBlueYellowRedScheme);

	//These are for internal use of the GUI program
	IMG::generic2DfieldPtr_t dynamicData; //to mediate between simulation and renderer data
	IMG::floats2Dfield_t floatField = IMG::createFloats2Dfield(width, height);
	dynamicData.storeFloatsField(&floatField);

	//TODO: next version of fViz2D should make this unecessary
	COLOR::rgbaF_t clear = COLOR::CLEAR;

	LOG_DEBUG("Starting renderer...");
	F_V2::rendererRetCode_st retCode = F_V2::rendererRetCode_st::STILL_RUNNING;
	std::thread renderThread = F_V2::spawnRendererOnNewThread(&dynamicData, &retCode, &clear, 
		                                                      &menuList, filenameFunc,
		                                                      &scheme, "Phase Field Model: CH", 
		                                                      768, 704, true, F_V2::defaultSynchCallback,
		                                                      PFM_GUI::pfmGuiBannerPathFromBinary);
	F_V2::saveCurrentImage(); //save initial state

	//This is the internal loop to keep the renderer data mirror fed
	float* data_ptr = floatField.data.get(); //the data in dynamicData as an array of floats
	size_t elements = floatField.size.getTotalElements();
	double lastReadTotalAbsoluteChangeSinceSaving = 0;
	while (retCode == F_V2::rendererRetCode_st::STILL_RUNNING) {
		//TODO: there should be better ways to do this : )
		for (size_t i = 0; i < elements; i++) {
			data_ptr[i] = (float)PFM::getActiveFieldConstPtr()->getElement(i);
		}

		if (checks_ptr->totalAbsoluteChangeSinceLastSave < lastReadTotalAbsoluteChangeSinceSaving) {
			//the total change was reset: a save was triggered: let's also save the renderer's image
			F_V2::saveCurrentImage();
		}
		lastReadTotalAbsoluteChangeSinceSaving = checks_ptr->totalAbsoluteChangeSinceLastSave;

		AZ::hybridBusySleepForMicros(std::chrono::microseconds(500));
	}

	renderThread.join();
	LOG_DEBUG("Render thread joined");

	//From here on we're winding down, closing stuff and giving some diagnostics:
	LOG_TRACE("Simulation being topped. Will join its thread");

	int stepsRan = PFM::stopSimulation();
	printf("Ran for %d steps\n", stepsRan);

	LOG_TRACE("Will save");

	//TODO: saving data on exit should be optional
	PFM::saveFieldData(true, true, true);

	LOG_INFO("Results saved");

	//Figure the return code and return
	bool result = (retCode == F_V2::rendererRetCode_st::OK);
	if (result) { LOG_INFO("OK"); }
	else { LOG_ERROR("Renderer or user-reported error"); }

	GETCHAR_FORCE_PAUSE;

	return result;
}