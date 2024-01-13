#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/miscDefines.hpp"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "fViz2D/API/FV2_API.hpp"
#include "fViz2D/API/GUI_API.hpp"

#include "PFM_API.hpp"
#include "PFM_defaults.hpp"
#include "PFM_tests.hpp"

#include "guiTests.hpp"
#include "GUIcontrolMain.hpp"
#include "guiMenus.hpp"

#define RUN_GUI_TESTS (-1)
#define SAVE_PARAMETERS_BEFORE_EACH_CHECK 1
#define OVERRIDE_DEFAULT_AND_TRY_TO_START_PAUSED 1 //unless specified by command line argument
#define SKIP_PAUSE_WHEN_SAVING_GUI_IMAGE 0

double g_changePerElementPerStepToStop = PFM::defaultAbsChangePerStepToStop;
uint64_t g_maximumSteps = PFM::defaulMaxSteps;
int g_checksAtChangeTreshold = 0;

namespace PFM_GUI { 
	bool g_restartSimulationAfterStopped = false;
	bool g_dtLoweredForFirstSteps = false;
    double g_originalDt = -1.0;
}

//TODO: A LOT OF THE PARSING COULD BE HANDLED VIA THE SIMULATION PROJECT

bool isArgumentDefault(int argument, char **argv) {
	return (strcmp(PFM::defaultArgument, argv[argument]) == 0);
}

void printArgumentsList() {
	for (int i = 0; i < PFM::mainsArgumentList::TOTAL_ARGS; i++) {
		printf("%s ", PFM::argumentNames[i]);
	}
	puts("\n");
}

//TODO: Sanitize input : )
//TODO: Deal with repetition between this and CLcontrolMain.cpp
bool processClInput(int* simToRun_ptr, PFM::simParameters_t* params_ptr, 
	                PFM::simConfig_t* config_ptr, int argc, char **argv) {

	LOG_TRACE("Will preccess command line inputs...");

	if (argc > 1 && argc != PFM::mainsArgumentList::TOTAL_ARGS) {
		LOG_ERROR("BAD NUBER OF CL ARGUMENTS: either just run the program or specify all arguments");
		printArgumentsList();
		return false;
	}

	*simToRun_ptr = (int)PFM::defaultSimToRun;
	//If arguments where provided, first figure out what simulation will be run:
	if(argc > 1) { 
		int index = (int)PFM::mainsArgumentList::SIM_TO_RUN;
		if (strcmp(PFM::defaultArgument, argv[index]) != 0) {
			//the defaultArgument wasn't passed for the simulation to run, so:
			*simToRun_ptr = atoi(argv[index]); 
		}
	}
		
	if(*simToRun_ptr == RUN_GUI_TESTS) { return true; } //nothing else is needed
	if ( *simToRun_ptr < 0 || *simToRun_ptr >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS ) { //bad sim id, quit
		LOG_ERROR("Simulation id requested through argument is not supported");
		return false;
	}	

	PFM::setSavingOnDATofTheParamsBeforeEachCheck(SAVE_PARAMETERS_BEFORE_EACH_CHECK);

	//Otherwise, lets first load the defaults for the chosen simulation
	*params_ptr = PFM::defaultSimParams[*simToRun_ptr];
	*config_ptr = PFM::defaultConfigs[*simToRun_ptr];

	if(OVERRIDE_DEFAULT_AND_TRY_TO_START_PAUSED) { config_ptr->startPaused = true; }

	//And then proceed to change any non-default values passed:
	for (int i = 1; i < argc; i++) {
		switch (i) {
		
			default: 
				LOG_ERROR("A valid argument is not being handled in processClInput");
				printf("arg: %d\n", i);
				return false;
		
			case PFM::mainsArgumentList::PROG_CALL: {
			} break;

			case PFM::mainsArgumentList::SIM_TO_RUN: {
			} break;

			case PFM::mainsArgumentList::LAMBDA: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->lambda));
			} break;

			case PFM::mainsArgumentList::GAMMA: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->gamma));
			} break;

			case PFM::mainsArgumentList::DT: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->dt));
			} break;

			case PFM::mainsArgumentList::CELLS: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				config_ptr->cells = atoi(argv[i]);
			} break;

			case PFM::mainsArgumentList::WIDTH: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				config_ptr->width = atoi(argv[i]);
			} break;

			case PFM::mainsArgumentList::HEIGHT: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				config_ptr->height = atoi(argv[i]);
			} break;

			case PFM::mainsArgumentList::INITIAL_COND: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				int cond = atoi(argv[i]);
				if (cond < 0 || cond >= (int)PFM::initialConditions::TOTAL_INITIAL_CONDS) {
					LOG_ERROR("Bad initial Condition");
					return false;
				}
				config_ptr->initialContidion = (PFM::initialConditions)cond;
			} break;

			case PFM::mainsArgumentList::BIAS: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(config_ptr->bias));
			} break;

			case PFM::mainsArgumentList::SEED: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				config_ptr->initialSeed = strtoll(argv[i], NULL, 10);
			} break;

			case PFM::mainsArgumentList::METHOD: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				int method = atoi(argv[i]);
				if (method < 0 || method >= (int)PFM::integrationMethods::TOTAL_METHODS) {
					LOG_ERROR("Bad Integration Method");
					return false;
				}
				config_ptr->method = (PFM::integrationMethods)method;
			} break;

			case PFM::mainsArgumentList::START_PAUSED: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				config_ptr->startPaused = atoi(argv[i]);
			} break;
			
			case PFM::mainsArgumentList::CHANGE_PER_ELEMENT_PER_STEP_TO_STOP: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &g_changePerElementPerStepToStop);
			} break;

			case PFM::mainsArgumentList::MAXIMUM_STEPS: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				g_maximumSteps = strtoll(argv[i], NULL, 10);
			} break;

			case PFM::mainsArgumentList::STEPS_PER_CHECK: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				int64_t steps = strtol(argv[i], NULL, 10);
				bool valid = (steps > 0) && (steps <= UINT32_MAX);
				if(!valid) { LOG_ERROR("Bad number of steps per check"); return false;}

				PFM::setMaxStepsPerCheckAdded((uint32_t)steps);
			} break;

			case PFM::mainsArgumentList::ABSOLUTE_CHANGE_PER_CHECK: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				double change = sscanf(argv[i], "%lf", &g_changePerElementPerStepToStop);
				if(change <= 0) { LOG_ERROR("Bad change per check"); return false; }

				PFM::setMaxTotalChangePerElementPerCheckAdded(change);
			} break;
				
			case PFM::mainsArgumentList::CALLER_KEY: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				PFM::setCallerKey((uint32_t)strtoll(argv[i], NULL, 10));
			} break;		

			case PFM::mainsArgumentList::ADAPTATIVE_DT: {
				if (strcmp(PFM::defaultArgument, argv[i]) == 0) { break; }
				params_ptr->useAdaptativeDt = atoi(argv[i]);
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
	
	do {
		PFM_GUI::g_restartSimulationAfterStopped = false;
		PFM_GUI::g_dtLoweredForFirstSteps = false;
		PFM_GUI::g_originalDt = params.dt;

		switch ((int)simToRun) {
			case RUN_GUI_TESTS:
				result = PFM_GUI_TESTS::guiLinkingAndDependencyTests();
			break;

			case (int)PFM::simFuncEnum::SINGLE_LAYER_CH_SIM:
				result = PFM_GUI::runSimulationWithGUI(&params, &config, PFM::getDirAndFileName); 
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
	} while(PFM_GUI::g_restartSimulationAfterStopped);

	if(result) { LOG_INFO("All ok"); }
	else { LOG_ERROR("Errors found"); }

	GETCHAR_FORCE_PAUSE;
	return !result; //0 for ok
}

//Returns empty menuList in case something goes wrong
void populateMenuList(GUI::menuDefinitionList_t* menuList_ptr, PFM::simConfig_t* simConfigFromMain_ptr, 
		                                                   PFM::simParameters_t* simParamsFromMain_ptr) {

	auto checks_ptr = PFM::getCheckDataPtr();
	auto config_ptr = PFM::getSimConfigPtr();
	auto params_ptr = PFM::getSimParamsPtr();

	if (checks_ptr == nullptr || config_ptr == nullptr || params_ptr == nullptr) {
		return;
	}

	menuList_ptr->push_back(PFM_GUI::getChecksMenuDefinition(checks_ptr));
	menuList_ptr->push_back(PFM_GUI::getConfigAndParamsMenuDefinition(config_ptr, params_ptr));
	menuList_ptr->push_back(PFM_GUI::getcontrolFlowMenuDefinition(simConfigFromMain_ptr, simParamsFromMain_ptr));

	return;
}

//Should also receive the number of steps to run and/or minimumAbsoluteChange
//Should receive stepsPerSave and/or totalAbsoluteChangePerSave as well
//Make sure if a positive number of steps was passed, manually closing the renderer still stop the simulation
//Add option for data saves to also save the image from the GUI
bool PFM_GUI::runSimulationWithGUI(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
	                               GUI::filenameCallback_func* filenameFunc, int stepsToRun, bool saveOnExit) {
	
	if ((int)config_ptr->simulFunc >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS) {
		LOG_ERROR("Invalid simulation function selected");
		return false;
	}

	if(config_ptr->simulFunc == PFM::simFuncEnum::DATA_CONTROL_TEST) { 
		LOG_WARN("The menu to validate this test is not supported on the current version of the controller. Will say that the manual inspection failed. TODO: Reimplement support");
	}

	LOG_DEBUG("Will run the simulation and display on the GUI");

	//Set up and run the simulation:
	PFM::setIntermediateBINsaves(true);
	PFM::initializeSimulation(*config_ptr);
	LOG_INFO("Simulation initialized. Will run the simulation");
	PFM::runForSteps(stepsToRun, *parameters_ptr, *config_ptr);
	
	//Form the menus:
	GUI::menuDefinitionList_t menuList;
	populateMenuList(&menuList, config_ptr, parameters_ptr);
	if (menuList.size() == 0) {
		LOG_ERROR("Could not form the menus");
		return false;
	}	

	//Loading the color scheme:
	//This should also be received
	COLOR::colorInterpolation_t scheme;
	scheme.loadScheme(&COLOR::defaultBlueYellowRedScheme);

	//These are for internal use of the GUI program
	IMG::generic2DfieldPtr_t dynamicData; //to mediate between simulation and renderer data
	IMG::floats2Dfield_t floatField = IMG::createFloats2Dfield(config_ptr->width, config_ptr->height);
	dynamicData.storeFloatsField(&floatField);

	//TODO: next version of fViz2D should make this unecessary
	COLOR::rgbaF_t clear = COLOR::CLEAR;

	LOG_DEBUG("Starting renderer...");
	F_V2::rendererRetCode_st retCode = F_V2::rendererRetCode_st::STILL_RUNNING;
	std::thread renderThread = F_V2::spawnRendererOnNewThread(&dynamicData, &retCode, &clear, 
		                                                      &menuList, filenameFunc,
		                                                      &scheme, "Phase Field Model: CH", 
		                                                      1000, 704, true, F_V2::defaultSynchCallback,
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

		if (PFM::getCheckDataPtr()->totalAbsoluteChangeSinceLastSave < lastReadTotalAbsoluteChangeSinceSaving) {
			//the total change was reset: a save was triggered: let's also save the renderer's image
			if(!SKIP_PAUSE_WHEN_SAVING_GUI_IMAGE) { PFM::pauseSimulation(); }
			F_V2::saveCurrentImage();
			if(!SKIP_PAUSE_WHEN_SAVING_GUI_IMAGE) { PFM::resumeSimulation(); }
		}
		lastReadTotalAbsoluteChangeSinceSaving = PFM::getCheckDataPtr()->totalAbsoluteChangeSinceLastSave;

		AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI/2));
	}

	renderThread.join();
	LOG_DEBUG("Render thread joined");

	//From here on we're winding down, closing stuff and giving some diagnostics:
	LOG_TRACE("Simulation being topped. Will join its thread");

	int stepsRan = PFM::stopSimulation();
	printf("Ran for %d steps\n", stepsRan);

	if(saveOnExit) { 
		LOG_TRACE("Will save");
		PFM::saveFieldData(true, true, true); 
		LOG_INFO("Results saved");
	}

	//Figure the return code and return
	bool result = (retCode == F_V2::rendererRetCode_st::OK);
	if (result) { LOG_INFO("OK"); }
	else { LOG_ERROR("Renderer or user-reported error"); }

	GETCHAR_FORCE_PAUSE;

	return result;
}