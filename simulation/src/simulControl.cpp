#include "PFM_API.hpp"
#include "PFM_tests.hpp"
#include "simulControl.hpp"

#include "fAux/API/timeHelpers.hpp"
#include "fAux/API/miscDefines.hpp"

using namespace PFM;

static SimulationControl controller;

///Definitions of the methods of the claas SimulationControl:

bool PFM::SimulationControl::isInitialized() const {
	return m_hasInitialized;
}

PeriodicDoublesLattice2D* PFM::SimulationControl::getFieldPtr() const {
	return m_activeLattice_ptr.get();
}
        
void PFM::SimulationControl::releaseField() {
	if(m_activeLattice_ptr != NULL) { m_activeLattice_ptr.release(); }
}

void PFM::SimulationControl::reinitializeController(fieldDimensions_t dimensions, uint32_t numberCells, 
	                                                                              double cellSeedValue) {
	
	if(isSimulationRunning()) { stop(); }
	releaseField();

	std::vector<double> initialData;
	size_t elements = dimensions.totalElements();
	initialData.reserve(elements);

	//The data to be used on initialization:
	int cellSpacing = elements;
	if(numberCells > 1) { cellSpacing /= numberCells; }
	for (size_t i = 0; i < elements; i++) {
		//To "seed" numberCells equally spaced cells with initialValue, and all others with zero:
		initialData.push_back(cellSeedValue * ( (i/cellSpacing) == (i/(double)cellSpacing) ) );
	}

	m_activeLattice_ptr = 
		std::unique_ptr<PeriodicDoublesLattice2D>(new PeriodicDoublesLattice2D(dimensions, initialData));

	m_cells = numberCells;
	m_shouldStop = false;
	m_stepsRan = 0;
    m_stepsToRun = 0;

	m_hasInitialized = true;
}

void PFM::SimulationControl::nonBlockingStop() {
	m_shouldStop = true;
}

void PFM::SimulationControl::stepsEnded() {
	m_isRunning = false;
	m_shouldStop = false;
}

int PFM::SimulationControl::stop() {
	nonBlockingStop();
	while(m_isRunning) {
		AZ::hybridBusySleepForMicros(std::chrono::microseconds(MS_TO_WAIT * MICROS_IN_A_MILLI));
	}

	if(m_stepsThread.joinable()) { m_stepsThread.join(); }

	m_isRunning = false;
	m_shouldStop = false;

	return m_stepsRan;
}

//TODO: an actual reasonable save system : p
bool PFM::SimulationControl::saveFieldToFile() const {
	if (!m_hasInitialized || !m_activeLattice_ptr->isInitialized() || !m_activeLattice_ptr->hasAllocated()) {
		return false;
	}

	auto dimensions = m_activeLattice_ptr->getFieldDimensions();

	std::string baseFilename = "sim" + std::to_string((int)m_lastSimulFuncUsed) + "_" + std::to_string(m_cells) 
		                    + "_" + std::to_string(dimensions.width) + "_" + std::to_string(dimensions.height) 
		                    + "_" + std::to_string(m_stepsRan) + "_" + std::to_string(m_initialSeed);
	
	FILE* fp_pgm = fopen((baseFilename + ".pgm").c_str(), "wb");
	if(fp_pgm == NULL) { return false; }

	FILE* fp_bin = fopen((baseFilename + ".bin").c_str(), "wb");
	if(fp_bin == NULL) { return false; }

	int maxColor = 255;
	fprintf(fp_pgm, "P5\n%d %d\n%d\n", (int)dimensions.width, (int)dimensions.height, maxColor);

	double value;
	for (int j = 0; j < (int)dimensions.height; j++) {
		for (int i = 0; i < (int)dimensions.width; i++) {
		  value = m_activeLattice_ptr->getDataPoint({i,j});
		  fprintf(fp_pgm, "%c", (char)(value*maxColor));
		  fprintf(fp_bin, "%f", value);
		}
	}
	fclose(fp_pgm);
	fclose(fp_bin);

	return true;
}

bool PFM::SimulationControl::checkIfShouldStop() {
	if(m_shouldStop || (m_stepsToRun > 0 && m_stepsRan >= m_stepsToRun)) {
		m_shouldStop = true;
	}

	return m_shouldStop;
}

int PFM::SimulationControl::getNumberCells() const {
	return m_cells;
}

double PFM::SimulationControl::getLastCellSeedValue() const {
	return m_lastCellSeedValue;
}

void PFM::SimulationControl::runForSteps(int steps, PFM::simFuncEnum simulationToRun) {
	if(controller.isSimulationRunning()) { return; }
	if((int)simulationToRun >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS) { return; }

	m_stepsToRun = steps;
	m_lastSimulFuncUsed = simulationToRun;
	m_isRunning = true;

	m_stepsThread = std::thread(*(simFunctionsPtrs_arr[(int)simulationToRun]), this, &m_stepsRan, &m_isRunning);
	m_stepsThread.detach();
}

bool PFM::SimulationControl::isSimulationRunning() const {
	return m_isRunning;
}

int PFM::SimulationControl::stepsAlreadyRan() const {
	return m_stepsRan;
}

void PFM::SimulationControl::resetStepsAlreadyRan() {
	m_stepsRan = 0;
}

///These API calls are really just wrappers to calls to methods of the controller:

PeriodicDoublesLattice2D* PFM::initializeSimulation(fieldDimensions_t dimensions, uint32_t numberCells) {
	
	controller.reinitializeController(dimensions, numberCells);
	return controller.getFieldPtr();
}

bool PFM::isSimulationRunning() {
	return controller.isSimulationRunning();
}

int PFM::stopSimulation() {
	return controller.stop();
}

int PFM::getStepsRan() {
	return controller.stepsAlreadyRan();
}

void PFM::resetStepsRan() {
	controller.resetStepsAlreadyRan();
}

bool PFM::saveFieldToFile() {
	return controller.saveFieldToFile();
}
	
void PFM::runForSteps(int stepsToRun, PFM::simFuncEnum simulationToRun) {
	controller.runForSteps(stepsToRun, simulationToRun);
}

PeriodicDoublesLattice2D* PFM::getActiveFieldPtr() {
	return controller.getFieldPtr();
}