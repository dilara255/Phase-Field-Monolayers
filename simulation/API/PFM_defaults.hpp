#pragma once

#include "PFM_data.hpp"
#include "PFM_API_enums.hpp"

namespace PFM {

	enum mainsArgumentList { PROG_CALL, SIM_TO_RUN, LAMBDA, GAMMA, DT, 
		                     CELLS, WIDTH, HEIGHT, INITIAL_COND, BIAS, SEED, METHOD, START_PAUSED,
							 CHANGE_PER_ELEMENT_PER_STEP_TO_STOP, MAXIMUM_STEPS,
						     STEPS_PER_CHECK, ABSOLUTE_CHANGE_PER_CHECK, CALLER_KEY,
	                         TOTAL_ARGS};
	static const char* argumentNames[TOTAL_ARGS] = { "ProgramCall", "SimFuncToRun(uint)", "Lamba(double)", "Gamma(double)",
		                                      "dt(double)", "Cells(uint)", "Width(uint)", "Height(uint)", 
											  "InitialCondition(uint)", "Bias(double)", "Seed(uint64)", 
											  "Method(uint)", "StartPaused(bool)", 
											  "ChangePerElementPerStepToStop(double)", "MaximumSteps(uint64)",
	                                          "StepsPerCheck(uint)", "ChangePerCheck(double)",
		                                      "CallerKey(uint)"};
	static const char* deafultArgument = "default";

	static const simFuncEnum defaultSimToRun = simFuncEnum::SINGLE_LAYER_CH_SIM;

	static const simConfig_t defaultConfigs[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {

		{0, 0, 0, 0, 0.0, 0, simFuncEnum::DATA_CONTROL_TEST, 
		 initialConditions::TOTAL_INITIAL_CONDS, 0.0, integrationMethods::TOTAL_METHODS, false, false},

		{0, 0, 128, 128, 1.0, DEFAULT_PRNG_SEED0, simFuncEnum::SINGLE_LAYER_CH_SIM, 
		 initialConditions::LINEAR_RANDOM, -0.3, integrationMethods::FTCS, false, false},

		{0, 50, 128, 128, 1.0, DEFAULT_PRNG_SEED0, simFuncEnum::MULTI_LAYER_CH_SIM, 
		 initialConditions::EVENLY_SPACED_INDEX, 0.0, integrationMethods::FTCS, true, false}
	};

	static const simParameters_t defaultSimParams[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {
		//{1.0, 7.824813, 0.043986, ...}
		{1.0, 3.0, 0.06, -1.0, -1.0}, {1.0, 3.0, 0.06, -1.0, -1.0}, {1.0, 3.0, 0.06, -1.0, -1.0}
	};

	static const uint64_t completelyArbitraryStepToUnlockFullDt = 50;
	//This is used to help avoid having way too many save in the first few frames
	static const double absoluteChangeRemainingFactor = 0.2;

	static const uint32_t defaulStepsPerCheck = 5000;
	static const double defaultAbsChangePerCheck = 0.0025;

	static const double defaultAbsChangePerStepToStop = 0.0000000015;
	static const uint64_t defaulMaxSteps = 6000000;
	static const int stepsAtChangeThresholdToStop = 100; 

	static const uint32_t defaultCallerKey = 0;
}