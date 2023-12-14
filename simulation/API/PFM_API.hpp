#pragma once

//TODO: make sure no threads can be orphaned under normal circunstances
//Also, thread safety is a crutch : )

#include "PFM_API_enums.hpp"
#include "PFM_data.hpp"
#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/prng.hpp"

namespace PFM {

	//If anything fails, returns NULL. Otherwise, returns a const pointer to the active field simulated
	PFM_API const PeriodicDoublesLattice2D* initializeSimulation(PFM::fieldDimensions_t dimensions, 
																	   uint32_t numberCells, 
																	   PFM::initialConditions initialCond = 
																			initialConditions::EVENLY_SPACED_INDEX,
		                                                               double bias = 0, 
																	   bool perCellLayer = false, 
	                                                                   uint64_t seed = DEFAULT_PRNG_SEED0);

	PFM_API bool isSimulationRunning();
	PFM_API int stopSimulation();
	PFM_API int getStepsRan();
	PFM_API void resetStepsRan();
	PFM_API bool saveFieldToFile();
	//Uses parameters and etc to build a default file-name. Compatible with GUI save button if sent as callback.
    PFM_API std::string getFileName(int steps, bool calledFromGui);
	
	//If steps <= 0, will run until manually stopped
	//If the simulation is already running or the simFuncEnum or method is bad, will do nothing (*no warning*)
	PFM_API void runForSteps(int stepsToRun, double lambda = 3.0, double gamma = 0.06, double dt = 1.0,
		                     PFM::simFuncEnum simulationToRun = PFM::simFuncEnum::SINGLE_LAYER_CH_SIM,
		                     integrationMethods method = integrationMethods::FTCS);
	
	//Returns NULL in case no field is active or the simulation hasn't been initialized
	PFM_API PeriodicDoublesLattice2D* getActiveFieldPtr();
	PFM_API const PeriodicDoublesLattice2D* getActiveFieldConstPtr();

	//These will return nullptr if the controller is not initialized
	PFM_API checkData_t* getCheckDataPtr();
	PFM_API const simConfig_t* getSimConfigPtr();
	PFM_API simParameters_t* getSimParamsPtr();

	//Update internal parameters bases on (potentially) new "external" parameters, sent previously by the user
	//Does nothing in case the controller is not initialized
	PFM_API void updatePhysicalParameters();
}