#pragma once

#include <memory>
#include <thread>
#include <functional>

#include "fAux/API/prng.hpp"

#include "PFM_API.hpp"
#include "PFM_data.hpp"

#include "fAux/API/prng.hpp"

#define CELL_SEED_VAL (1.0) //used to mark a point as a cell to be expanded. TODO: deprecate?
#define MS_TO_WAIT 10
#define DEFAULT_STEPS_PER_CHECK 500
#define DEFAULT_ABSOLUTE_CHANGE_PER_CHECK (0.01)

namespace PFM {

    //*******************SimulationFunctions*************************
    //These are the definitions of the possible simulation functions:
    class SimulationControl;
    typedef void SimulationSteps_fn(SimulationControl* controller_ptr, uint64_t* stepCount_ptr, 
                                    const bool* shouldPause_ptr, bool* isRunning_ptr, integrationMethods method);
    
    //WARNING: these should always be kept in synch with PFM::simFuncEnum (@PFM_API.hpp)
    SimulationSteps_fn dataAndControllerTest_fn;
    SimulationSteps_fn singleLayerCHsim_fn;
    SimulationSteps_fn multiLayerCHsim_fn;
    static SimulationSteps_fn* simFunctionsPtrs_arr[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {
        &dataAndControllerTest_fn, &singleLayerCHsim_fn, &multiLayerCHsim_fn
    };
    //***************************************************************


    //********************SimulationControl**************************
    //This class holds the data and references necessary to control the simulation
    class SimulationControl {

    public:

        bool isInitialized() const;
        CurrentAndLastPerioricDoublesLattice2D* getRotatingBaseFieldPtr();
        PeriodicDoublesLattice2D* getBaseFieldPtr();
        PeriodicDoublesLattice2D* getLastDphisAndTempKsFieldPtr();
        
        std::vector<std::unique_ptr<PeriodicDoublesLattice2D>>* getLayerFieldsVectorPtr() const;
    
        //Returns a pointer to the active base field
        //May be a pointer to a distinc base-field or to the last layer of the rotating field
        //By default, points to the base lattice
        PeriodicDoublesLattice2D* getActiveFieldPtr();
        //Which one is set by setRotatingLastAsActive() and setBaseAsActive()
        PeriodicDoublesLattice2D* setRotatingLastAsActive(); //also returns the new active pointer
        PeriodicDoublesLattice2D* setBaseAsActive(); //also returns the new active pointer

        checkData_t* getActiveFieldsCheckDataPtr();

        //If not yet initialized, initializes and creates a new field
        //Otherwise, destroys the old field, reinitializes and creates a new field with "dimensions"
        //Defaults to initialConditions::EVENLY_SPACED_INDEX if a bad condition is passed
        void reinitializeController(PFM::simConfig_t config);

        //Spawns a new thread which will run the simulation.
        //Thread is joined either when the steps are over or from a call to stop() / nonBlockingStop().
        //Does nothing in case the simulation is already running or a bad simFuncEnum is passed.
        void runForSteps(uint64_t steps,  double lambda, double gamma, double dt, simFuncEnum simulationToRun, 
                         integrationMethods methodToUse = integrationMethods::FTCS);
        //Asks the controller to stop and waits for confirmation (in MS_TO_WAIT ms sleep cycles)
        //When the simulation actually stops, its thread is joined. Does nothing if the simulation isn't running.
        //Returns the amount of steps ran
        int stop();

        void pause();
        void resume();
        bool shouldBePaused() const;

        //Copies the data from the currentStep of the rotating fields into the baseField
        void mirrorRotatingOnBase();
        //Copies the data from the baseField into the currentStep of the rotating fields
        void mirrorBaseOnRotating();

        //Returns false in case the field was unitialized or unallocated
        bool saveFieldData(bool savePGM, bool saveBIN, bool saveDAT) const;
        //Same as above, but uses the controllers values to decide what is saved or not
        //Set using the "setIntermediateXXXsaves" methods - they're all false by default (this does nothing then)
        bool saveFieldData() const;
        void setAused(double newA);
        void setKused(double newK);
        void setLambdaUsed(double newLambda);
        void setGammaUsed(double newGamma);
        //Just returns in case the new lambda = 0
        void updatePhysicalParametersFromInternals();
        void setDTused(double newK);

        void setMaxStepsPerCheckAdded(size_t newStepsPerCheckSaved);
        //Set to zero in case newMaxTotalChangePerCheck < 0 (will save every step)
        void setMaxTotalChangePerElementPerCheckAdded(double newMaxTotalChangePerCheck);

	    void setIntermediateDATsaves(bool shouldSave);
	    void setIntermediatePGMsaves(bool shouldSave);
	    void setIntermediateBINsaves(bool shouldSave);
        void setSavingOnDATofTheParamsBeforeEachCheck(bool shouldSave);

        const simConfig_t* getLastSimConfigPtr() const;
        std::string getSimDataString() const;

        simParameters_t* getLastSimParametersPtr();
        std::string getSimParamsString() const;

        size_t getActiveFieldsCheckVectorElements() const;
        //Returns an empty string in case the checkNumber is bad
        std::string getActiveFieldsChecksString(size_t checkNumber) const;
        //Returns an empty string in case the checkNumber is bad
        std::string getActiveFieldsParamStringBeforeAGivenCheck(size_t checkNumber) const;
        
        void printSimDataAndParams() const;

        bool isSimulationRunning() const;
        const bool* getIsPaused_ptr() const;
        bool checkIfShouldStop();
        int getNumberCells() const;
        double getLastCellSeedValue() const;
        bool shouldStillExpandSeeds() const;
        int getStepsPerCheckSaved() const;
        double getAbsoluteChangePerCheckSaved() const;

        uint64_t stepsAlreadyRan() const;
        void resetStepsAlreadyRan();

        void updateEpochTimeSimCall();

    private:

        void releaseFields();
        void stepsEnded();
        //Asks the controller to stop, but otherwise keeps going. *Use with caution*.
        //WARNING: When the simulation actually stops, its thread is STILL HANGING. Does nothing if not running.
        void nonBlockingStop();
        //These will just return in case their dividend is zero:
        void updateGammaLambda();
        void updateKandA();

        bool m_hasInitialized = false;
        bool m_isRunning = false;
        bool m_shouldStop = false;
        bool m_shouldBePaused = false;
        int m_stepsToRun = 0;
        bool m_seedsNeedExpanding = false;
        int m_stepsPerCheckSaved = DEFAULT_STEPS_PER_CHECK;
        double m_absoluteChangePerCheckSaved = DEFAULT_ABSOLUTE_CHANGE_PER_CHECK;

        bool m_saveDATonIntermediateChecks = false;
        bool m_savePGMonIntermediateChecks = false;
        bool m_saveBINonIntermediateChecks = false;
        bool m_saveOnDATtheParamsBeforeEachCheck = false;

        std::thread m_stepsThread;

        simConfig_t m_simConfigs;
        simParameters_t m_simParameters;
        
        //To be filled with actual data:
        std::unique_ptr<PeriodicDoublesLattice2D> m_baseLattice_ptr;
        std::unique_ptr<PFM::CurrentAndLastPerioricDoublesLattice2D> m_rotatingBaseLattice_ptr;
        std::unique_ptr<PeriodicDoublesLattice2D> m_lastDphisAndTempKsField_ptr;

        PeriodicDoublesLattice2D* m_activeBaseField_ptr;

        std::vector<std::unique_ptr<PeriodicDoublesLattice2D>> m_perCellLaticePtrs;
    };
    //***************************************************************
}