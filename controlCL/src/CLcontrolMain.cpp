#include "fAux/API/logAPI.hpp"
#include "fAux/API/timeHelpers.hpp"
#include "fAux/API/miscDefines.hpp"

#include "PFM_API.hpp"
#include "PFM_defaults.hpp"
#include "PFM_tests.hpp"
#include "CLcontrolMain.hpp"

#define RUN_CLI_TESTS -1                            
#define DEFAULT_CHANGE_PER_ELEMENT_PER_STEP_TO_STOP (0.0000000015)
#define DEFAULT_MAX_STEPS (6000000)
#define STEPS_AT_CHANGE_THRESHOLD_TO_ACTUALLY_STOP 1

double g_changePerElementPerStepToStop = DEFAULT_CHANGE_PER_ELEMENT_PER_STEP_TO_STOP;
uint64_t g_maximumSteps = DEFAULT_MAX_STEPS;
int g_checksAtChangeTreshold = 0;

bool isArgumentDefault(int argument, char **argv) {
	return (strcmp(PFM_CLI::deafultArgument, argv[argument]) == 0);
}

void printArgumentsList() {
	for (int i = 0; i < PFM_CLI::mainsArgumentList::TOTAL_ARGS; i++) {
		printf("%s ", PFM_CLI::argumentNames[i]);
	}
	puts("\n");
}

//TODO: Sanitize input : )
bool processClInput(int* simToRun_ptr, PFM::simParameters_t* params_ptr, 
	                PFM::simConfig_t* config_ptr, int argc, char **argv) {

	LOG_TRACE("Will preccess command line inputs...");

	if (argc > 1 && argc != PFM_CLI::mainsArgumentList::TOTAL_ARGS) {
		LOG_ERROR("BAD NUBER OF CL ARGUMENTS: either just run the program or specify all arguments");
		printArgumentsList();
		return false;
	}

	*simToRun_ptr = (int)PFM::defaultSimToRun;
	//If arguments where provided, first figure out what simulation will be run:
	if(argc > 1) { 
		int index = (int)PFM_CLI::mainsArgumentList::SIM_TO_RUN;
		if (strcmp(PFM_CLI::deafultArgument, argv[index]) != 0) {
			//the defaultArgument wasn't passed for the simulation to run, so:
			*simToRun_ptr = atoi(argv[index]); 
		}
	}
		
	if(*simToRun_ptr == RUN_CLI_TESTS) { return true; } //nothing else is needed
	if ( *simToRun_ptr < 0 || *simToRun_ptr >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS ) { //bad sim id, quit
		LOG_ERROR("Simulation id requested through argument is not supported");
		return false;
	}	

	//Otherwise, lets first load the defaults for the chosen simulation
	*params_ptr = PFM::defaultSimParams[*simToRun_ptr];
	*config_ptr = PFM::defaultConfigs[*simToRun_ptr];

	//And then proceed to change any non-default values passed:
	for (int i = 1; i < argc; i++) {
		switch (i) {
		
			default: 
				LOG_ERROR("A valid argument is not being handled in processClInput");
				printf("arg: %d\n", i);
				return false;
		
			case PFM_CLI::mainsArgumentList::PROG_CALL: {
			} break;

			case PFM_CLI::mainsArgumentList::SIM_TO_RUN: {
			} break;

			case PFM_CLI::mainsArgumentList::LAMBDA: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->lambda));
			} break;

			case PFM_CLI::mainsArgumentList::GAMMA: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->gamma));
			} break;

			case PFM_CLI::mainsArgumentList::DT: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->dt));
			} break;

			case PFM_CLI::mainsArgumentList::CELLS: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->cells = atoi(argv[i]);
			} break;

			case PFM_CLI::mainsArgumentList::WIDTH: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->width = atoi(argv[i]);
			} break;

			case PFM_CLI::mainsArgumentList::HEIGHT: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->height = atoi(argv[i]);
			} break;

			case PFM_CLI::mainsArgumentList::INITIAL_COND: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				int cond = atoi(argv[i]);
				if (cond < 0 || cond >= (int)PFM::initialConditions::TOTAL_INITIAL_CONDS) {
					LOG_ERROR("Bad initial Condition");
					return false;
				}
				config_ptr->initialContidion = (PFM::initialConditions)cond;
			} break;

			case PFM_CLI::mainsArgumentList::BIAS: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(config_ptr->bias));
			} break;

			case PFM_CLI::mainsArgumentList::SEED: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->initialSeed = strtoll(argv[i], NULL, 10);
			} break;

			case PFM_CLI::mainsArgumentList::METHOD: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				int method = atoi(argv[i]);
				if (method < 0 || method >= (int)PFM::integrationMethods::TOTAL_METHODS) {
					LOG_ERROR("Bad Integration Method");
					return false;
				}
				config_ptr->method = (PFM::integrationMethods)method;
			} break;

			case PFM_CLI::mainsArgumentList::START_PAUSED: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->startPaused = atoi(argv[i]);
			} break;

			case PFM_CLI::mainsArgumentList::CHANGE_PER_ELEMENT_PER_STEP_TO_STOP: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &g_changePerElementPerStepToStop);
			} break;

			case PFM_CLI::mainsArgumentList::MAXIMUM_STEPS: {
				if (strcmp(PFM_CLI::deafultArgument, argv[i]) == 0) { break; }
				g_maximumSteps = strtoll(argv[i], NULL, 10);
			} break;
		}
	}

	return true;
}

int main(int argc, char **argv) {

	PFM::simConfig_t config;
	PFM::simParameters_t params;
	int simToRun;

	if(!processClInput(&simToRun, &params, &config, argc, argv)) { 
		LOG_ERROR("Bad input. Quitting");
		GETCHAR_FORCE_PAUSE;
		return 1; 
	}

	bool result = false;
	LOG_INFO("Input processed. Will call the selected simulation function...");
	printf("maxSteps = %llu\nchangePerElementPerStepToStop: %f\n", g_maximumSteps, g_changePerElementPerStepToStop);

	switch ((int)simToRun) {
		case RUN_CLI_TESTS:
			result = PFM::linkingTest();
		break;

		case (int)PFM::simFuncEnum::SINGLE_LAYER_CH_SIM:
			result = PFM_CLI::runSimulationFromCL(&params, &config); 
		break;

		case (int)PFM::simFuncEnum::MULTI_LAYER_CH_SIM:
			LOG_WARN("MULTI LAYER SIM TEMPORARILY DISABLED. Todo: re-enable");
		break;

		case (int)PFM::simFuncEnum::DATA_CONTROL_TEST:
			LOG_WARN("DATA CONTROL TESTS TEMPORARILY DISABLED. Todo: re-enable");
		break;

		default:
			LOG_ERROR("Chosen simulation not supported by the GUI program");
		break;
	}

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

bool shouldStop(PFM::checkData_t* check_ptr) {
	if ((check_ptr->lastAbsoluteChange / check_ptr->stepsDuringLastCheckPeriod) <= g_changePerElementPerStepToStop) {
		g_checksAtChangeTreshold++;
	}
	return (check_ptr->stepsAtLastCheck >= g_maximumSteps) || 
		   (g_checksAtChangeTreshold >= STEPS_AT_CHANGE_THRESHOLD_TO_ACTUALLY_STOP);
}

//TODO: Pass stuff in just like to the GUI counterpart
bool PFM_CLI::runSimulationFromCL(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
		                                                             int stepsToRun, bool saveOnExit) {
	
	if ((int)config_ptr->simulFunc >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS) {
		LOG_ERROR("Invalid simulation function selected");
		return false;
	}

	if(config_ptr->simulFunc == PFM::simFuncEnum::DATA_CONTROL_TEST) { 
		LOG_WARN("The menu to validate this test is not supported on the current version of the controller. Will say that the manual inspection failed. TODO: Reimplement support");
	}

	LOG_DEBUG("Will run the simulation from the CLI");

	//Set up and run the simulation:
	PFM::setIntermediateBINsaves(true);
	PFM::initializeSimulation(*config_ptr);
	LOG_INFO("Simulation initialized. Will run the simulation");
	PFM::runForSteps(stepsToRun, *parameters_ptr, *config_ptr);
	
	auto check_ptr = PFM::getCheckDataPtr();
	g_checksAtChangeTreshold = 0; //in case restart or something like taht is implemented
	uint64_t stepsRan = PFM::getStepsRan(); //"initial" step
	while (PFM::isSimulationRunning()) {
		//Wait for a new step in the simulation:
		while(stepsRan == PFM::getStepsRan()) { 
			AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI/2));
		}
		//Update current step:
		stepsRan = PFM::getStepsRan(); 
		//And check if we should stop:
		if(shouldStop(check_ptr)) { LOG_INFO("Stop criteria reached"); stepsRan = PFM::stopSimulation(); }
	}

	LOG_DEBUG("Simulation ended and thread joined");
	printf("Ran for %d steps\n", stepsRan);

	if(saveOnExit) { 
		LOG_TRACE("Will save");
		PFM::saveFieldData(true, true, true); 
		LOG_INFO("Results saved");
	}

	GETCHAR_FORCE_PAUSE;

	return true;
}