#include "PFM_API.hpp"
#include "fViz2D/prototypes.hpp"

namespace PFM_GUI {

	//stepsToRun <= 0: run until manually stopped
	bool runSimulationWithGUI(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
		                      GUI::filenameCallback_func* filenameFunc, int stepsToRun = -1, bool saveOnExit = true);
}