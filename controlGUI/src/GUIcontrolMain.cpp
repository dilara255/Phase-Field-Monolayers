#include "fViz2D/API/FV2_API.hpp"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"

#include "guiTests.hpp"
#include "guiControlMain.hpp"

#define RUN_GUI_TESTS 0
#define RUN_SIM_DATACTRL_TESTS 0
#define RUN_SINGLE_LAYER_CH_SIM_NO_OLD 0
#define RUN_SINGLE_LAYER_CH_SIM_WITH_OLD 1
#define RUN_MULTI_LAYER_CH_SIM 0
#define TOTAL_SIM_FUNCS ((int)PFM::simFuncEnum::TOTAL_SIM_FUNCS)

typedef struct parameters_st {
	int width, height, cells;
	PFM::initialConditions initialCond = PFM::initialConditions::EVENLY_SPACED_INDEX;
} parameters_t;

parameters_t defaultParamsPerSimulType[TOTAL_SIM_FUNCS] = {
	{512, 512, 50},
	{128, 128, 1, PFM::initialConditions::BALANCED_RANDOM},
	{128, 128, 1, PFM::initialConditions::BALANCED_RANDOM},
	{128, 128, 5}
};

int main() {

	bool result = true;

	if(RUN_GUI_TESTS) { result &= PFM_GUI_TESTS::guiLinkingAndDependencyTests(); }
	if(RUN_SIM_DATACTRL_TESTS) { result &= PFM_GUI::runSimulation(PFM::simFuncEnum::DATA_CONTROL_TEST); }
	if(RUN_SINGLE_LAYER_CH_SIM_NO_OLD) { 
		result &= PFM_GUI::runSimulation(PFM::simFuncEnum::SINGLE_LAYER_CH_SIM_NO_OLD); 
	}
	if(RUN_SINGLE_LAYER_CH_SIM_WITH_OLD) { 
		result &= PFM_GUI::runSimulation(PFM::simFuncEnum::SINGLE_LAYER_CH_SIM_WITH_OLD); 
	}
	if(RUN_MULTI_LAYER_CH_SIM) { result &= PFM_GUI::runSimulation(PFM::simFuncEnum::MULTI_LAYER_CH_SIM); }

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

bool PFM_GUI::runSimulation(PFM::simFuncEnum simulationFunctionToRun) {
	if ((int)simulationFunctionToRun >= TOTAL_SIM_FUNCS) {
		LOG_ERROR("Invalid simulation function selected");
		return false;
	}

	LOG_DEBUG("Will run the simulation and display on the GUI");

	int simIndex = (int)simulationFunctionToRun;

	int width = defaultParamsPerSimulType[simIndex].width;
	int height = defaultParamsPerSimulType[simIndex].height;
	int cells = defaultParamsPerSimulType[simIndex].cells;

	PFM::fieldDimensions_t dimensions = {(size_t)width, (size_t)height};

	
	auto dataField_ptr = PFM::initializeSimulation(dimensions, cells, 
		                                           defaultParamsPerSimulType[simIndex].initialCond);

	LOG_INFO("Simulation initialized");
		
	IMG::floats2Dfield_t floatField = IMG::createFloats2Dfield(width, height);
	
	F_V2::rendererRetCode_st retCode = F_V2::rendererRetCode_st::STILL_RUNNING;

	//These are just for compatibility with the test version of fViz2D
	//TODO: update fViz2D version and get rid of this : )
	COLOR::rgbaF_t clear = COLOR::CLEAR;
	COLOR::rgbaF_t tint = COLOR::FULL_WHITE;
	//****************************************************************

	bool works = false;

	LOG_DEBUG("Starting renderer...");
	std::thread renderThread = F_V2::spawnRendererOnNewThreadF(&works, &floatField, &clear, &tint, &retCode);

	LOG_INFO("Will run the simulation");
	PFM::runForSteps(-1, simulationFunctionToRun);

	size_t elements = floatField.size.getTotalElements();

	while (retCode == F_V2::rendererRetCode_st::STILL_RUNNING) {
		for (size_t i = 0; i < elements; i++) {
			floatField.data[i] = (float)dataField_ptr->getElement(i);
		}

		AZ::hybridBusySleepForMicros(std::chrono::microseconds(1000));
	}

	LOG_INFO("Simulation ended");
		
	LOG_TRACE("Will stop the thread");

	renderThread.join();
	LOG_DEBUG("Render thread joined");

	int stepsRan = PFM::stopSimulation();

	LOG_DEBUG("Simulation thread joined");
	printf("Ran for %d steps\n", stepsRan);

	LOG_TRACE("Will save");

	PFM::saveFieldToFile();

	LOG_INFO("Results saved");

	bool result = (retCode == F_V2::rendererRetCode_st::OK && works);
	if (result) { LOG_INFO("OK"); }
	else { LOG_ERROR("Renderer or user-reported error"); }

	GETCHAR_FORCE_PAUSE;

	return result;
}