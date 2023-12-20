#pragma once

#include "PFM_data.hpp"
#include "PFM_API_enums.hpp"

//DATA_CONTROL_TEST, SINGLE_LAYER_CH_SIM, MULTI_LAYER_CH_SIM

namespace PFM {

	static const simFuncEnum defaultSimToRun = simFuncEnum::SINGLE_LAYER_CH_SIM;

	static const simConfig_t defaultConfigs[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {

		{0, 0, 0, 0, 0.0, 0, simFuncEnum::DATA_CONTROL_TEST, 
		 initialConditions::TOTAL_INITIAL_CONDS, 0.0, integrationMethods::TOTAL_METHODS, false},

		{0, 0, 128, 128, 1.0, DEFAULT_PRNG_SEED0, simFuncEnum::SINGLE_LAYER_CH_SIM, 
		 initialConditions::LINEAR_RANDOM, -0.3, integrationMethods::FTCS, false},

		{0, 50, 128, 128, 1.0, DEFAULT_PRNG_SEED0, simFuncEnum::MULTI_LAYER_CH_SIM, 
		 initialConditions::EVENLY_SPACED_INDEX, 0.0, integrationMethods::FTCS, true}
	};

	static const simParameters_t defaultSimParams[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {
		//{1.0, 7.824813, 0.043986, ...}
		{1.0, 3.0, 0.06, -1.0, -1.0}, {1.0, 3.0, 0.06, -1.0, -1.0}, {1.0, 3.0, 0.06, -1.0, -1.0}
	};
}