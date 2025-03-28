#pragma once

namespace PFM {

	enum class simFuncEnum { DATA_CONTROL_TEST, SINGLE_LAYER_CH_SIM, MULTI_LAYER_CH_SIM, TOTAL_SIM_FUNCS};
	//TODO: Add BEST_CIRCULAR_PACKING, BEST_HEXAGONAL_PACKING and FORCED_CONCENTRATION_RANDOM
	enum class initialConditions { EVENLY_SPACED_INDEX, LINEAR_RANDOM, TOTAL_INITIAL_CONDS};
	enum class integrationMethods { FTCS, FTCS_WITH_SUBS, HEUN, RK2, RK4, TOTAL_METHODS };

	enum mainsArgumentList { PROG_CALL, SIM_TO_RUN, LAMBDA, GAMMA, DT, 
		                    CELLS, WIDTH, HEIGHT, INITIAL_COND, BIAS, SEED, METHOD, START_PAUSED,
							CHANGE_PER_ELEMENT_PER_STEP_TO_STOP, MAXIMUM_STEPS,
						    STEPS_PER_CHECK, ABSOLUTE_CHANGE_PER_CHECK, CALLER_KEY, ADAPTATIVE_DT,
							MAX_CHANGE_PER_STEP, MAX_SPEEDUP_MULT, MIN_SLOWDOWN_MULT, SAFE_MAX_DT,
	                        TOTAL_MAIN_ARGS};
}