#include "PFM_API.hpp"
#include "fViz2D/prototypes.hpp"

namespace PFM_GUI {

	enum mainsArgumentList { PROG_CALL, SIM_TO_RUN, LAMBDA, GAMMA, DT, 
		                     CELLS, WIDTH, HEIGHT, INITIAL_COND, BIAS, SEED, METHOD, START_PAUSED,
							 CHANGE_PER_ELEMENT_PER_STEP_TO_STOP, MAXIMUM_STEPS,
						     STEPS_PER_CHECK, ABSOLUTE_CHANGE_PER_CHECK,
	                         TOTAL_ARGS};
	const char* argumentNames[TOTAL_ARGS] = { "ProgramCall", "SimFuncToRun(uint)", "Lamba(double)", "Gamma(double)",
		                                      "dt(double)", "Cells(uint)", "Width(uint)", "Height(uint)", 
											  "InitialCondition(uint)", "Bias(double)", "Seed(uint64)", 
											  "Method(uint)", "StartPaused(bool)", 
											  "ChangePerElementPerStepToStop(double)", "MaximumSteps(uint64)",
	                                          "StepsPerCheck(uint)", "ChangePerCheck(double)"};
	const char* deafultArgument = "default";

	//stepsToRun <= 0: run until manually stopped
	bool runSimulationWithGUI(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
		                      GUI::filenameCallback_func* filenameFunc, int stepsToRun = -1, bool saveOnExit = true);
}