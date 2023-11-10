#include "fAux/API/logAPI.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"
#include "CLcontrolMain.hpp"

#define RUN_CLI_TESTS 0

int main() {
	if(RUN_CLI_TESTS) { PFM::linkingTest();	return 0;}
	else { return !PFM_CLI::runSimulation(); }
}

bool PFM_CLI::runSimulation() {

	LOG_DEBUG("Will run the simulation from the CLI and save the results as an image");

	int width = 512;
	int height = 512;
	int cells = 50;
	int stepsToRun = 166;

	PFM::fieldDimensions_t dimensions = {(size_t)width, (size_t)height};

	PFM::initializeSimulation(dimensions, cells);

	LOG_INFO("Simulation initialized");
		
	PFM::runForSteps(stepsToRun);

	LOG_TRACE("Simulation running");

	while (PFM::isSimulationRunning()) {
		AZ::hybridBusySleepForMicros(std::chrono::microseconds(1000));
	}

	LOG_TRACE("Simulation ended. Will stop the thread");

	int stepsRan = PFM::stopSimulation();
	printf("Ran for %d steps.\n", stepsRan);

	LOG_TRACE("Will save");

	PFM::saveFieldToFile();

	LOG_INFO("Results saved");

	bool result = (stepsRan == stepsToRun);
	if (result) { LOG_INFO("OK"); }
	else { LOG_ERROR("A different number of steps was expected"); }

	GETCHAR_FORCE_PAUSE;

	return !result;

}