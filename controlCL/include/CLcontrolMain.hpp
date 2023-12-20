#include "PFM_API.hpp"

namespace PFM_CLI {

	enum mainsArgumentList { PROG_CALL, SIM_TO_RUN, LAMBDA, GAMMA, DT, 
		                     CELLS, WIDTH, HEIGHT, INITIAL_COND, BIAS, SEED, METHOD, START_PAUSED,
							 CHANGE_PER_ELEMENT_PER_STEP_TO_STOP, MAXIMUM_STEPS,
	                         TOTAL_ARGS};
	const char* argumentNames[TOTAL_ARGS] = { "ProgramCall", "SimFuncToRun(uint)", "Lamba(double)", "Gamma(double)",
		                                      "dt(double)", "Cells(uint)", "Width(uint)", "Height(uint)", 
											  "InitialCondition(uint)", "Bias(double)", "Seed(uint64)", 
											  "Method(uint)", "StartPaused(bool)", 
											  "ChangePerElementPerStepToStop(double)", "MaximumSteps(uint64)"};

	const char* deafultArgument = "default";

	bool runSimulationFromCL(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
		                                            int stepsToRun = -1, bool saveOnExit = true);
}