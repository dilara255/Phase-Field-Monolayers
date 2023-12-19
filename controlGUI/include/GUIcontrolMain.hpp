#include "PFM_API.hpp"
#include "fViz2D/prototypes.hpp"

namespace PFM_GUI {

	enum mainsArgumentList { PROG_CALL, SIM_TO_RUN, LAMBDA, GAMMA, DT, 
		                     CELLS, WIDTH, HEIGHT, INITIAL_COND, BIAS, SEED, METHOD,
	                         TOTAL_ARGS};
	const char* argumentNames[TOTAL_ARGS] = { "ProgramCall", "SimFuncToRun(uint)", "Lamba(double)", "Gamma(double)",
		                                      "dt(double)", "Cells(uint)", "Width(uint)", "Height(uint)", 
											  "InitialCondition(uint)", "Bias(double)", "Seed(uint64)", "Method(uint)"};
	const char* deafultArgument = "default";

	//stepsToRun <= 0: run until manually stopped
	bool runSimulationWithGUI(PFM::simParameters_t parameters, PFM::simConfig_t config,
		                      GUI::filenameCallback_func* filenameFunc, int stepsToRun = -1, bool saveOnExit = true);
}