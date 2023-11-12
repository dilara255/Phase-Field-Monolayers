#pragma once

#include <memory>
#include <thread>
#include <functional>

#include "fAux/API/prng.hpp"

#include "PFM_API.hpp"
#include "PFM_data.hpp"

#define CELL_SEED_VAL (1.0)
#define MS_TO_WAIT 10

namespace PFM {

    //*******************SimulationFunctions*************************
    //These are the definitions of the possible simulation functions:
    class SimulationControl;
    typedef void SimulationSteps_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr);
    
    //WARNING: these should always be kept in synch with PFM::simFuncEnum (@PFM_API.hpp)
    SimulationSteps_fn dataAndControllerTest_fn;
    SimulationSteps_fn singleLayerSim_fn;
    static SimulationSteps_fn* simFunctionsPtrs_arr[(int)simFuncEnum::TOTAL_SIM_FUNCS] = {
        &dataAndControllerTest_fn, &singleLayerSim_fn
    };
    //***************************************************************


    //********************SimulationControl**************************
    //This class holds the data and references necessary to control the simulation
    class SimulationControl {
    //TODO: probably should be a singleton : )

    public:
        bool isInitialized() const;
        PeriodicDoublesLattice2D* getFieldPtr() const; //TODO: does this really work with const in the front?
        
        //If not yet initialized, initializes and creates a new field
        //Otherwise, destroys the old field, reinitializes and creates a new field with "dimensions"
        void reinitializeController(fieldDimensions_t dimensions, uint32_t numberCells, 
                                    bool perCellLayer, double cellSeedValue = CELL_SEED_VAL);

        //Spawns a new thread which will run the simulation.
        //Thread is joined either when the steps are over or from a call to stop() / nonBlockingStop().
        //Does nothing in case the simulation is already running or a bad simFuncEnum is passed.
        void runForSteps(int steps, simFuncEnum simulationToRun);
        //Asks the controller to stop and waits for confirmation (in MS_TO_WAIT ms sleep cycles)
        //When the simulation actually stops, its thread is joined. Does nothing if the simulation isn't running.
        //Returns the amount of steps ran
        int stop();

        //Returns false in case the field was unitialized or unallocated
        bool saveFieldToFile() const;
        void setAused(double newA);
        void setKused(double newK);
        void setDTused(double newK);

        bool isSimulationRunning() const;
        bool checkIfShouldStop();
        int getNumberCells() const;
        double getLastCellSeedValue() const;

        int stepsAlreadyRan() const;
        void resetStepsAlreadyRan();

    private:
            
        void releaseField();
        void stepsEnded();
        //Asks the controller to stop, but otherwise keeps going. *Use with caution*.
        //WARNING: When the simulation actually stops, its thread is STILL HANGING. Does nothing not running.
        void nonBlockingStop();

        bool m_hasInitialized = false;
        bool m_isRunning = false;
        bool m_shouldStop = false;
        int m_stepsRan = 0;
        int m_stepsToRun = 0;
        
        std::thread m_stepsThread;

        int m_cells = 0;
        double m_lastCellSeedValue = CELL_SEED_VAL;
        uint64_t m_initialSeed = DEFAULT_PRNG_SEED0;
        PFM::simFuncEnum m_lastSimulFuncUsed = PFM::simFuncEnum::TOTAL_SIM_FUNCS;        
        double m_lastK = -1;
        double m_lastA = -1;
        double m_lastDT = -1;
        
        std::unique_ptr<PeriodicDoublesLattice2D> m_baseLattice_ptr = NULL;
        std::vector<std::unique_ptr<PeriodicDoublesLattice2D>> m_perCellLaticePtrs;
    };
    //***************************************************************
}