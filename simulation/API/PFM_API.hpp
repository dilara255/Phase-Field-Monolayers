#pragma once

//TODO: make sure no threads can be orphaned under normal circunstances
//Also, thread safety is a crutch : )

#include "PFM_API_enums.hpp"
#include "PFM_data.hpp"
#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/prng.hpp"

namespace PFM {

	//If anything fails, returns NULL. Otherwise, returns a const pointer to the active field simulated
	PFM_API const PeriodicDoublesLattice2D* initializeSimulation(simConfig_t config);

	//If steps <= 0, will run until manually stopped
	//If the simulation is already running or the simFuncEnum or method is bad, will do nothing (*no warning*)
	PFM_API void runForSteps(int stepsToRun, simParameters_t parameters, simConfig_t config);
	PFM_API int stopSimulation();

	PFM_API bool isSimulationRunning();
	PFM_API int getStepsRan();
	PFM_API void resetStepsRan();

	PFM_API bool saveFieldData(bool savePGM, bool saveBIN, bool saveDAT);
	//Uses parameters and etc to build a default file-name. Compatible with GUI save button if sent as callback.
    PFM_API std::string getFileName(int steps, bool calledFromGui);

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

	//The values set by the following two functions are used to decide when to add a new checkData entry
	//The entry is triggered when either condition is met
	//Whenever the "DAT" data is saved, it includes a list of all the checkData entries so far
	//This one is compared to the number of steps since the last entry
	//If seto to zero, will add a new entry every step (not reccomended)
	PFM_API void setMaxStepsPerCheckAdded(size_t newMaxStepsPerCheck);
	//And this one is compared to the absolute change per element since the last entry (with spillover)
	//Will be set to 0 if a smaller value is passed: entries will be added every step (not reccomended)
	PFM_API void setMaxTotalChangePerElementPerCheckAdded(double newMaxTotalChangePerCheck);

	//These change wether intermediate saves of each data are made when each new checkData is added
	PFM_API void setIntermediateDATsaves(bool shouldSave); //The DATs have config, param and CUMMULATIVE check data
	PFM_API void setIntermediatePGMsaves(bool shouldSave); //The PGMs are simple 0-255 images of the network
	PFM_API void setIntermediateBINsaves(bool shouldSave); //The BINs are the binary data of the network
	//Default false. If true, each checkData entry will be preceeded by the parameters at the last check
	PFM_API void setSavingOnDATofTheParamsBeforeEachCheck(bool shouldSave);
}