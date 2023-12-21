#include "PFM_API.hpp"

namespace PFM_CLI {

	bool runSimulationFromCL(PFM::simParameters_t* parameters_ptr, PFM::simConfig_t* config_ptr,
		                                            int stepsToRun = -1, bool saveOnExit = true);
}