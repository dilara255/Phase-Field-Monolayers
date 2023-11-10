#include "fViz2D/API/FV2_API.hpp"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"

#include "guiTests.hpp"
#include "guiControlMain.hpp"

#define RUN_GUI_TESTS 0
#define RUN_SIM_DATACTRL_TESTS 0
#define RUN_SINGLE_LAYER_SIM 1

int main() {

	bool result = true;

	if(RUN_GUI_TESTS) { result &= PFM_GUI_TESTS::guiLinkingAndDependencyTests(); }
	if(RUN_SIM_DATACTRL_TESTS) { result &= PFM_GUI::runSimulation(PFM::simFuncEnum::DATA_CONTROL_TEST); }
	if(RUN_SINGLE_LAYER_SIM) { result &= PFM_GUI::runSimulation(PFM::simFuncEnum::SINGLE_LAYER_SIM); }

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

bool PFM_GUI::runSimulation(PFM::simFuncEnum simulationFunctionToRun) {

	LOG_DEBUG("Will run the simulation and display on the GUI");

	int width = 512;
	int height = 512;
	int cells = 50;

	PFM::fieldDimensions_t dimensions = {(size_t)width, (size_t)height};

	auto dataField_ptr = PFM::initializeSimulation(dimensions, cells);
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