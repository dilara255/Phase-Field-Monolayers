#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "fViz2D/API/FV2_API.hpp"
#include "fViz2D/API/GUI_API.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"

#include "guiTests.hpp"
#include "guiControlMain.hpp"
#include "guiMenus.hpp"

#define RUN_GUI_TESTS (-1)
#define TOTAL_SIM_FUNCS ((int)PFM::simFuncEnum::TOTAL_SIM_FUNCS)

//TODO: A LOT OF THE PARSING COULD BE HANDLED VIA THE SIMULATION PROJECT

bool isArgumentDefault(int argument, char **argv) {
	return (strcmp(PFM_GUI::deafultArgument, argv[argument]) == 0);
}

void printArgumentsList() {
	for (int i = 0; i < PFM_GUI::mainsArgumentList::TOTAL_ARGS; i++) {
		printf("%s ", PFM_GUI::argumentNames[i]);
	}
	puts("\n");
}

//TODO: Sanitize input : )
bool processClInput(int* simToRun_ptr, PFM::simParameters_t* params_ptr, 
	                             PFM::simConfig_t* config_ptr, int argc, char **argv) {

	LOG_TRACE("Will preccess command line inputs...");

	if (argc > 1 && argc != PFM_GUI::mainsArgumentList::TOTAL_ARGS) {
		LOG_ERROR("BAD NUBER OF CL ARGUMENTS: either just run the program or specify all arguments");
		printArgumentsList();
		return false;
	}

	//TODO: first figure out the simulation func, then load *its* defaults, then proceed
	params_ptr->loadDefaults();
	config_ptr->loadDefaults();
	*simToRun_ptr = (int)config_ptr->simulFunc;

	for (int i = 1; i < argc; i++) {
		switch (i) {
		
			default: 
				LOG_ERROR("Absurd error parsing command line input");
				return false;
		
			//Otherwise, the argument is treated here:
			case PFM_GUI::mainsArgumentList::SIM_TO_RUN: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				*simToRun_ptr = atoi(argv[i]);		
			} break;

			case PFM_GUI::mainsArgumentList::LAMBDA: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->lambda));
			} break;

			case PFM_GUI::mainsArgumentList::GAMMA: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->gamma));
			} break;

			case PFM_GUI::mainsArgumentList::DT: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				sscanf(argv[i], "%lf", &(params_ptr->dt));
			} break;

			case PFM_GUI::mainsArgumentList::CELLS: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->cells = atoi(argv[i]);
			} break;

			case PFM_GUI::mainsArgumentList::WIDTH: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->width = atoi(argv[i]);
			} break;

			case PFM_GUI::mainsArgumentList::HEIGHT: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->height = atoi(argv[i]);
			} break;

			case PFM_GUI::mainsArgumentList::INITIAL_COND: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				int cond = atoi(argv[i]);
				if (cond < 0 || cond >= (int)PFM::initialConditions::TOTAL_INITIAL_CONDS) {
					LOG_ERROR("Bad initial Condition");
					return false;
				}
				config_ptr->initialContidion = (PFM::initialConditions)cond;
			} break;

			case PFM_GUI::mainsArgumentList::BIAS: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				scanf(argv[i], "%lf", &(config_ptr->bias));
			} break;

			case PFM_GUI::mainsArgumentList::SEED: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				config_ptr->initialSeed = strtoll(argv[i], NULL, 10);
			} break;

			case PFM_GUI::mainsArgumentList::METHOD: {
				if (strcmp(PFM_GUI::deafultArgument, argv[i]) == 0) { break; }
				int method = atoi(argv[i]);
				if (method < 0 || method >= (int)PFM::integrationMethods::TOTAL_METHODS) {
					LOG_ERROR("Bad Integration Method");
					return false;
				}
				config_ptr->method = (PFM::integrationMethods)method;
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
	
	switch ((int)simToRun) {
		case RUN_GUI_TESTS:
			result = PFM_GUI_TESTS::guiLinkingAndDependencyTests();
		break;

		case (int)PFM::simFuncEnum::SINGLE_LAYER_CH_SIM:
			result = PFM_GUI::runSimulationWithGUI(params, config, PFM::getFileName); 
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

//Returns empty menuList in case something goes wrong
void populateMenuList(GUI::menuDefinitionList_t* menuList_ptr) {

	auto checks_ptr = PFM::getCheckDataPtr();
	auto config_ptr = PFM::getSimConfigPtr();
	auto params_ptr = PFM::getSimParamsPtr();

	if (checks_ptr == nullptr || config_ptr == nullptr || params_ptr == nullptr) {
		return;
	}

	menuList_ptr->push_back(PFM_GUI::getChecksMenuDefinition(checks_ptr));
	menuList_ptr->push_back(PFM_GUI::getConfigAndParamsMenuDefinition(config_ptr, params_ptr));

	return;
}

//TODO: receive config and param structs
//Should also receive the number of steps to run and/or minimumAbsoluteChange
//Should receive stepsPerSave and/or totalAbsoluteChangePerSave as well
//Make sure if a posite number of steps was passed, manually closing the renderer still stops the simulation
//Add option for data saves to also save the image from the GUI
//See config and params, receive those plus whatever
bool PFM_GUI::runSimulationWithGUI(PFM::simParameters_t parameters, PFM::simConfig_t config,
	                               GUI::filenameCallback_func* filenameFunc, int stepsToRun, bool saveOnExit) {
	if ((int)config.simulFunc >= TOTAL_SIM_FUNCS) {
		LOG_ERROR("Invalid simulation function selected");
		return false;
	}

	if(config.simulFunc == PFM::simFuncEnum::DATA_CONTROL_TEST) { 
		LOG_WARN("The menu to validate this test is not supported on the current version of the controller. Will say that the manual inspection failed. TODO: Reimplement support");
	}

	LOG_DEBUG("Will run the simulation and display on the GUI");

	//Set up and run the simulation:
	PFM::initializeSimulation(config);
	LOG_INFO("Simulation initialized. Will run the simulation");
	PFM::runForSteps(stepsToRun, parameters, config);
	
	//Form the menus:
	GUI::menuDefinitionList_t menuList;
	populateMenuList(&menuList);
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
	IMG::floats2Dfield_t floatField = IMG::createFloats2Dfield(config.width, config.height);
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

		if (PFM::getCheckDataPtr()->totalAbsoluteChangeSinceLastSave < lastReadTotalAbsoluteChangeSinceSaving) {
			//the total change was reset: a save was triggered: let's also save the renderer's image
			F_V2::saveCurrentImage();
		}
		lastReadTotalAbsoluteChangeSinceSaving = PFM::getCheckDataPtr()->totalAbsoluteChangeSinceLastSave;

		AZ::hybridBusySleepForMicros(std::chrono::microseconds(500));
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