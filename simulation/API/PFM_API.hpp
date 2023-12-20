#pragma once

//TODO: make sure no threads can be orphaned under normal circunstances
//Also, thread safety is a crutch : )

#include "PFM_API_enums.hpp"
#include "PFM_data.hpp"
#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/prng.hpp"

namespace PFM {

	//If anything fails, returns NULL. Otherwise, returns a const pointer to the active field simulated
	//In case the simulation is already running, will return a pointr to the current active field and do nothing
	//(even if paused - make sure to call stopSimulation() first, or check with isSimulationRunning());
	PFM_API const PeriodicDoublesLattice2D* initializeSimulation(simConfig_t config);

	//If steps <= 0, will run until manually stopped
	//If the simulation is already running or the simFuncEnum or method is bad, will do nothing (*no warning*)
	//Note that changes in config won't be applied to new calls to runForSteps unless initializeSimulation is called
	PFM_API void runForSteps(uint64_t stepsToRun, simParameters_t parameters, simConfig_t config);
	//In case the simulation is paused, it will first be resumed before actually stopping
	PFM_API int stopSimulation();
	PFM_API void pauseSimulation();
	PFM_API void resumeSimulation();

	//Returns true even if the simulation is actually paused
	PFM_API bool isSimulationRunning();
	//Only returns true if the simulation isSimulationRunning() is also true
	PFM_API bool isSimulationPaused();
	PFM_API uint64_t getStepsRan();
	PFM_API void resetStepsRan();

	//Returns false in case the field was unitialized or unallocated
	PFM_API bool saveFieldData(bool savePGM, bool saveBIN, bool saveDAT);
	//Same as above, but uses the controllers values to decide what is saved or not
    //Set using the "setIntermediateXXXsaves" methods - they're all false by default (this does nothing then)
	PFM_API bool saveFieldDataAccordingToController();
	//Uses parameters and etc to build a default file-name. Compatible with GUI save button if sent as callback.
	//If the directory doesnt exist, it is created. Does not include extension.
    PFM_API std::string getDirAndFileName(int steps, bool calledFromGui);

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
	PFM_API double getLambdaFromKandA(double k, double A);
	PFM_API double getGammaFromKandA(double k, double A);

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