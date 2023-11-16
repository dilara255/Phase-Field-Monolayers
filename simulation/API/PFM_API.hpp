#pragma once

//TODO: make sure no threads can be orphaned under normal circunstances
//Also, thread safety is a crutch : )

#include "PFM_data.hpp"

namespace PFM {

	enum class simFuncEnum { DATA_CONTROL_TEST, SINGLE_LAYER_CH_SIM_NO_OLD, SINGLE_LAYER_CH_SIM_WITH_OLD,
		                     MULTI_LAYER_CH_SIM, TOTAL_SIM_FUNCS};
	//TODO: Add BEST_CIRCULAR_PACKING, BEST_HEXAGONAL_PACKING e FORCED_CONCENTRATION_RANDOM
	enum class initialConditions { EVENLY_SPACED_INDEX, BALANCED_RANDOM, TOTAL_INITIAL_CONDS};

	//If anything fails, returns NULL. Otherwise, returns a const pointer to the active field simulated
	PFM_API const PeriodicDoublesLattice2D* initializeSimulation(PFM::fieldDimensions_t dimensions, 
																	   uint32_t numberCells, 
																	   PFM::initialConditions initialCond = 
																			initialConditions::EVENLY_SPACED_INDEX,
																	   bool perCellLayer = false);

	PFM_API bool isSimulationRunning();
	PFM_API int stopSimulation();
	PFM_API int getStepsRan();
	PFM_API void resetStepsRan();
	PFM_API bool saveFieldToFile();
	
	//If <= 0, will run until manually stopped
	//If the simulation is already running or the simFuncEnum is bad, will do nothing (*no warning*)
	PFM_API void runForSteps(int stepsToRun, PFM::simFuncEnum simulationToRun);
	
	//Returns NULL in case no field is active or the simulation hasn't been initialized
	PFM_API PeriodicDoublesLattice2D* getActiveFieldPtr();
	PFM_API const PeriodicDoublesLattice2D* getActiveFieldConstPtr();
}