#include "fAux/API/logAPI.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"
#include "CLcontrolMain.hpp"

#define RUN_CLI_TESTS 0
#define RUN_SIM_DATACTRL_TESTS 0
#define RUN_SINGLE_LAYER_CH_SIM 0
#define RUN_MULTI_LAYER_CH_SIM 0
#define TOTAL_SIM_FUNCS ((int)PFM::simFuncEnum::TOTAL_SIM_FUNCS)

int main() {
	bool result = true;

	if(RUN_CLI_TESTS) { PFM::linkingTest(); }
	if(RUN_SIM_DATACTRL_TESTS) { result &= PFM_CLI::runSimulationFromCL(PFM::simFuncEnum::DATA_CONTROL_TEST); }
	if(RUN_SINGLE_LAYER_CH_SIM) { result &= PFM_CLI::runSimulationFromCL(PFM::simFuncEnum::SINGLE_LAYER_CH_SIM); }
	if(RUN_MULTI_LAYER_CH_SIM) { result &= PFM_CLI::runSimulationFromCL(PFM::simFuncEnum::MULTI_LAYER_CH_SIM); }

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

//TODO: Pass stuff in just like to the GUI counterpart
bool PFM_CLI::runSimulationFromCL(PFM::simFuncEnum simulationFunctionToRun) {

	//TODO: change to reflect new saving
	LOG_DEBUG("Will run the simulation from the CLI and save the results as an image");

	//TODO: should check that there is a stop condition (either stepsToRun > 0 or minimumAbsoluteChange > 0)

	//vvvv This should be "standardized" with the GUi counterpart through SIM API vvvv
	int width = 512;
	int height = 512;
	int cells = 50;
	int stepsToRun = 166;

	PFM::fieldDimensions_t dimensions = {(size_t)width, (size_t)height};

	PFM::initializeSimulation(dimensions, cells);

	LOG_INFO("Simulation initialized");
		
	PFM::runForSteps(stepsToRun, 3.0, 0.06, 1, simulationFunctionToRun);

	LOG_TRACE("Simulation running");
	//^^^^ This should be "standardized" with the GUi counterpart through SIM API ^^^^

	while (PFM::isSimulationRunning()) {
		AZ::hybridBusySleepForMicros(std::chrono::microseconds(1000));
	}

	//From here on we're winding down, closing stuff and giving some diagnostics:
	LOG_TRACE("Simulation being topped. Will join its thread");

	int stepsRan = PFM::stopSimulation();
	printf("Ran for %d steps.\n", stepsRan);

	LOG_TRACE("Will save");

	//TODO: saving data on exit should be optional
	PFM::saveFieldData(true, true, true);

	LOG_INFO("Results saved");

	//Figure the return code and return
	bool result = (stepsRan == stepsToRun || stepsToRun <= 0);
	if (result) { LOG_INFO("OK"); }
	else { LOG_ERROR("A different number of steps was expected"); }

	GETCHAR_FORCE_PAUSE;

	return result;

}