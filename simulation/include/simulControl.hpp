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

namespace PFM {

    //*******************SimulationFunctions*************************
    //These are the definitions of the possible simulation functions:
    class SimulationControl;
    typedef void SimulationSteps_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr,
                                                                                     integrationMethods method);
    
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

        //TODO:
        //add an activeField_ptr. This will be the one the API will give the caller
        //add method to set getBaseFieldPtr to return pointer to m_rotatingBaseLattice_ptr's currentLayer
        //add method to set getBaseFieldPtr to return pointer to it's own field
        //These should also affect what the activeField_ptr actually points to

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
        void reinitializeController(fieldDimensions_t dimensions, uint32_t numberCells, 
                                    PFM::initialConditions initialCond, bool perCellLayer, 
                                    double bias = 0, double cellSeedValue = CELL_SEED_VAL,
                                    uint64_t prngSeed = DEFAULT_PRNG_SEED0);

        //Spawns a new thread which will run the simulation.
        //Thread is joined either when the steps are over or from a call to stop() / nonBlockingStop().
        //Does nothing in case the simulation is already running or a bad simFuncEnum is passed.
        void runForSteps(int steps, simFuncEnum simulationToRun, 
                         integrationMethods methodToUse = integrationMethods::FTCS);
        //Asks the controller to stop and waits for confirmation (in MS_TO_WAIT ms sleep cycles)
        //When the simulation actually stops, its thread is joined. Does nothing if the simulation isn't running.
        //Returns the amount of steps ran
        int stop();

        //Copies the data from the currentStep of the rotating fields into the baseField
        void mirrorRotatingOnBase();
        //Copies the data from the baseField into the currentStep of the rotating fields
        void mirrorBaseOnRotating();

        //Returns false in case the field was unitialized or unallocated
        bool saveFieldToFile() const;
        void setAused(double newA);
        void setKused(double newK);
        void setDTused(double newK);
        void setStepsPerCheckSaved(int newStepsPerCheckSaved);

        const simData_t* getLastSimDataPtr() const;
        std::string getSimDataString() const;

        const simParameters_t* getLastSimParametersPtr() const;
        std::string getSimParamsString() const;

        void printSimDataAndParams() const;

        bool isSimulationRunning() const;
        bool checkIfShouldStop();
        int getNumberCells() const;
        double getLastCellSeedValue() const;
        bool shouldStillExpandSeeds() const;
        int getStepsPerCheckSaved() const;

        int stepsAlreadyRan() const;
        void resetStepsAlreadyRan();

    private:
            
        void releaseFields();
        void stepsEnded();
        //Asks the controller to stop, but otherwise keeps going. *Use with caution*.
        //WARNING: When the simulation actually stops, its thread is STILL HANGING. Does nothing not running.
        void nonBlockingStop();
        void updateGammaLambda();

        bool m_hasInitialized = false;
        bool m_isRunning = false;
        bool m_shouldStop = false;
        int m_stepsToRun = 0;
        bool m_seedsNeedExpanding = false;
        int m_stepsPerCheckSaved = DEFAULT_STEPS_PER_CHECK;

        std::thread m_stepsThread;

        simData_t m_lastSimData;
        simParameters_t m_lastSimParameters;
        
        //To be filled with actual data:
        std::unique_ptr<PeriodicDoublesLattice2D> m_baseLattice_ptr;
        std::unique_ptr<PFM::CurrentAndLastPerioricDoublesLattice2D> m_rotatingBaseLattice_ptr;
        std::unique_ptr<PeriodicDoublesLattice2D> m_lastDphisAndTempKsField_ptr;

        PeriodicDoublesLattice2D* m_activeBaseField_ptr;

        std::vector<std::unique_ptr<PeriodicDoublesLattice2D>> m_perCellLaticePtrs;
    };
    //***************************************************************
}