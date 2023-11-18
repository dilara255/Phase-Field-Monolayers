#include "PFM_API.hpp"
#include "PFM_tests.hpp"
#include "simulControl.hpp"

#include "fAux/API/timeHelpers.hpp"
#include "fAux/API/miscDefines.hpp"
#include "fAux/API/prng.hpp"

using namespace PFM;

static SimulationControl controller;

///Definitions of the methods of the claas SimulationControl:

bool PFM::SimulationControl::isInitialized() const {
	return m_hasInitialized;
}

CurrentAndLastPerioricDoublesLattice2D* PFM::SimulationControl::getRotatingBaseFieldPtr() {
	return m_rotatingBaseLattice_ptr.get();
}

PeriodicDoublesLattice2D* PFM::SimulationControl::getBaseFieldPtr() {
	return m_baseLattice_ptr.get();
}

PeriodicDoublesLattice2D* PFM::SimulationControl::getLastDphisAndTempKsFieldPtr() {
	return m_lastDphisAndTempKsField_ptr.get();
}

std::vector<std::unique_ptr<PeriodicDoublesLattice2D>>* PFM::SimulationControl::getLayerFieldsVectorPtr() const {
	return (std::vector<std::unique_ptr<PeriodicDoublesLattice2D>>*)&m_perCellLaticePtrs;
}

PeriodicDoublesLattice2D* PFM::SimulationControl::getActiveFieldPtr() {
	return m_activeBaseField_ptr;
}

PeriodicDoublesLattice2D* PFM::SimulationControl::setRotatingLastAsActive() {
	m_activeBaseField_ptr = m_rotatingBaseLattice_ptr->getPointerToLast();
	return m_activeBaseField_ptr;
}

PeriodicDoublesLattice2D* PFM::SimulationControl::setBaseAsActive() {
	m_activeBaseField_ptr = m_baseLattice_ptr.get();
	return m_activeBaseField_ptr;	
}

checkData_t* PFM::SimulationControl::getActiveFieldsCheckDataPtr() {
	return &m_activeBaseField_ptr->checks;
}
        
void PFM::SimulationControl::releaseFields() {
	
	m_activeBaseField_ptr = NULL;
	if (m_rotatingBaseLattice_ptr != NULL) {
		m_rotatingBaseLattice_ptr->releaseFields();
		m_rotatingBaseLattice_ptr.release();
	}
	if (m_baseLattice_ptr != NULL) {
		m_baseLattice_ptr.release();
	}
	if (m_lastDphisAndTempKsField_ptr != NULL) {
		m_lastDphisAndTempKsField_ptr.release();
	}

	m_hasInitialized = false;
}

//TODO: this whole thing should be reimplemented, and probably split apart a bit
void PFM::SimulationControl::reinitializeController(fieldDimensions_t dimensions, uint32_t numberCells, 
													PFM::initialConditions initialCond, bool perCellLayer, 
	                                                                                 double cellSeedValue) {
	
	if(isSimulationRunning()) { stop(); }
	releaseFields();

	size_t elements = dimensions.totalElements();
	
	m_lastDphisAndTempKsField_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
			new PeriodicDoublesLattice2D(dimensions, ALL_CELLS_ID)
	);

	std::vector<double> initialData;

	initialData.reserve(elements);

	//The data to be used on initialization:
	int cellSpacing = elements;
	if(numberCells > 1) { cellSpacing /= numberCells; }
	bool centerSingleCell = false; //TODO: pass as parameter

	if(!perCellLayer) {
		if (initialCond == PFM::initialConditions::EVENLY_SPACED_INDEX) {
			for (size_t i = 0; i < elements; i++) {
				//To "seed" numberCells equally spaced cells with initialValue, and all others with zero:
				if (numberCells == 1 && centerSingleCell) {
					initialData.push_back(cellSeedValue * ( i == (elements + dimensions.width)/2) );
				}
				else {
					initialData.push_back(cellSeedValue * ( (i/cellSpacing) == (i/(double)cellSpacing) ) );
				}
			}
			m_seedsNeedExpanding = true;
		}
		else {
			uint64_t seed = m_initialSeed;
			for (size_t i = 0; i < elements; i++) {
				//Each elemnt starts with a random value between -0.5 and 1.5:
				double value = -0.5 + 2*(AZ::draw1spcg32(&seed)/(double)UINT32_MAX);
				initialData.push_back(value);
			}
			m_seedsNeedExpanding = false;
		}
		
		m_rotatingBaseLattice_ptr = std::unique_ptr<CurrentAndLastPerioricDoublesLattice2D>(
			new CurrentAndLastPerioricDoublesLattice2D(dimensions, ALL_CELLS_ID, initialData)
		);
		m_baseLattice_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
			new PeriodicDoublesLattice2D(dimensions, ALL_CELLS_ID, initialData)
		);
	}
	else {
		//TODO: Implement other seedings for this
		m_rotatingBaseLattice_ptr = std::unique_ptr<CurrentAndLastPerioricDoublesLattice2D>(
			new CurrentAndLastPerioricDoublesLattice2D(dimensions, ALL_CELLS_ID)
		);
		m_baseLattice_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
			new PeriodicDoublesLattice2D(dimensions, ALL_CELLS_ID, initialData)
		);

		m_perCellLaticePtrs.reserve(numberCells);
		for (uint32_t cell = 0; cell < numberCells; cell++) {
			for (size_t i = 0; i < elements; i++) {
				//To "seed" numberCells equally spaced cells with initialValue, and all others with zero:
				initialData.push_back((double)(i == cellSpacing * cell));
			}

			m_perCellLaticePtrs.push_back(
				std::unique_ptr<PeriodicDoublesLattice2D>(
					new PeriodicDoublesLattice2D(dimensions, cell, initialData)
				)
			);
			initialData.clear();
		}
	}

	m_cells = numberCells;
	m_shouldStop = false;
	m_stepsRan = 0;
    m_stepsToRun = 0;
	setBaseAsActive();

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

void PFM::SimulationControl::mirrorRotatingOnBase() {
	m_baseLattice_ptr.get()->mirrorAllDataFrom(m_rotatingBaseLattice_ptr.get()->getPointerToCurrent());
}

void PFM::SimulationControl::mirrorBaseOnRotating() {
	m_rotatingBaseLattice_ptr.get()->getPointerToCurrent()->mirrorAllDataFrom(m_baseLattice_ptr.get());
}

//TODO: an actual reasonable save system : p
bool PFM::SimulationControl::saveFieldToFile() const {

	if (!m_hasInitialized || m_activeBaseField_ptr == NULL) { return false; }
	
	if(!m_activeBaseField_ptr->isInitialized() || !m_activeBaseField_ptr->hasAllocated()) {
		return false;
	}

	auto dimensions = m_activeBaseField_ptr->getFieldDimensions();

	std::string baseFilename = "sim" + std::to_string((int)m_lastSimulFuncUsed) 
		                     + "_A" + std::to_string(m_lastA) + "_k" + std::to_string(m_lastK)
		                     + "_dt" + std::to_string(m_lastDT)
		                     + "_" + std::to_string(m_cells) + "_" + std::to_string(dimensions.width) 
		                     + "_" + std::to_string(dimensions.height) + "_" + std::to_string(m_stepsRan) 
		                     + "_" + std::to_string(m_initialSeed);
	
	FILE* fp_pgm = fopen((baseFilename + ".pgm").c_str(), "wb");
	if(fp_pgm == NULL) { return false; }

	FILE* fp_bin = fopen((baseFilename + ".bin").c_str(), "wb");
	if(fp_bin == NULL) { return false; }

	int maxColor = 255;
	fprintf(fp_pgm, "P5\n%d %d\n%d\n", (int)dimensions.width, (int)dimensions.height, maxColor);

	double value;
	for (int j = 0; j < (int)dimensions.height; j++) {
		for (int i = 0; i < (int)dimensions.width; i++) {
		  value = m_activeBaseField_ptr->getDataPoint({i,j});
		  fprintf(fp_pgm, "%c", (char)(value*maxColor));
		  fprintf(fp_bin, "%f", value);
		}
	}
	fclose(fp_pgm);
	fclose(fp_bin);

	return true;
}

void PFM::SimulationControl::setAused(double newA) {
	m_lastA = newA;
}

void PFM::SimulationControl::setKused(double newK) {
	m_lastK = newK;
}

void PFM::SimulationControl::setDTused(double newDT) {
	m_lastDT = newDT;
}

void PFM::SimulationControl::setStepsPerCheckSaved(int newStepsPerCheckSaved) {
	m_stepsPerCheckSaved = newStepsPerCheckSaved;
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

bool PFM::SimulationControl::shouldStillExpandSeeds() const {
	return m_seedsNeedExpanding;
}

int PFM::SimulationControl::getStepsPerCheckSaved() const {
	return m_stepsPerCheckSaved;
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

const PeriodicDoublesLattice2D* PFM::initializeSimulation(fieldDimensions_t dimensions, 
	                                                      uint32_t numberCells, 
	                                                      PFM::initialConditions initialCond, 
													      bool perCellLayer) {
	
	controller.reinitializeController(dimensions, numberCells, initialCond, perCellLayer);
	return controller.getActiveFieldPtr();
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
	return controller.getActiveFieldPtr();
}

const PeriodicDoublesLattice2D* PFM::getActiveFieldConstPtr() {
	return (const PeriodicDoublesLattice2D*)controller.getActiveFieldPtr();
}