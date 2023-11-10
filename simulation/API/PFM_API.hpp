#pragma once

//TODO: make sure no threads can be orphaned under normal circunstances
//Also, thread safety is a crutch : )

#include "PFM_data.hpp"

namespace PFM {

	//If anything fails, returns NULL. Otherwise, returns the field to be used by the simulation
	PFM_API PeriodicDoublesLattice2D* initializeSimulation(PFM::fieldDimensions_t dimensions, 
		                                                                uint32_t numberCells);

	PFM_API bool isSimulationRunning();
	PFM_API int stopSimulation();
	PFM_API int getStepsRan();
	PFM_API void resetStepsRan();
	PFM_API bool saveFieldToFile();
	
	//If <= 0, will run until manually stopped
	PFM_API void runForSteps(int stepsToRun);

	//Returns NULL in case no field is active or the simulation hasn't been initialized
	PFM_API PeriodicDoublesLattice2D* getActiveFieldPtr();
}