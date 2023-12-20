#include <algorithm>

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
	                                                double bias, double cellSeedValue, uint64_t prngSeed) {
	
	if(isSimulationRunning()) { stop(); }
	releaseFields();

	m_lastSimData.initialSeed = prngSeed;
	m_lastSimData.width = dimensions.width;
	m_lastSimData.height = dimensions.height;
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
			uint64_t seed = m_lastSimData.initialSeed;
			for (size_t i = 0; i < elements; i++) {
				//Each elemnt starts with a random value between -0.5 and 1.5:
				double value = bias -0.5 + 2*(AZ::draw1spcg32(&seed)/(double)UINT32_MAX);
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

	m_lastSimData.cells = numberCells;
	m_shouldStop = false;
	m_lastSimData.stepsRan = 0;
    m_stepsToRun = 0;
	setBaseAsActive();

	m_lastSimData.initialContidion = initialCond;
	m_lastSimData.bias = bias;

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

	return m_lastSimData.stepsRan;
}

void PFM::SimulationControl::mirrorRotatingOnBase() {
	m_baseLattice_ptr.get()->mirrorAllDataFrom(m_rotatingBaseLattice_ptr.get()->getPointerToCurrent());
}

void PFM::SimulationControl::mirrorBaseOnRotating() {
	m_rotatingBaseLattice_ptr.get()->getPointerToCurrent()->mirrorAllDataFrom(m_baseLattice_ptr.get());
}

std::string PFM::getFileName(int steps, bool calledFromGUI) {
	
	//gather the relevant data:
	auto dimensions = controller.getActiveFieldPtr()->getFieldDimensions();
	const PFM::simConfig_t* simData_ptr = controller.getLastSimConfigPtr();
	const PFM::simParameters_t* simParams_ptr = controller.getLastSimParametersPtr();

	std::string fileName = "";	
	if(calledFromGUI) { fileName += "m"; }

	return  "Sim" + std::to_string((int)simData_ptr->simulFunc) 
		    + "m" + std::to_string((int)simData_ptr->method)
		    + "_ini" + std::to_string((int)simData_ptr->initialContidion) 
			+ "_b" + std::to_string(simData_ptr->bias)
		    + "_A" + std::to_string(simParams_ptr->A) + "_k" + std::to_string(simParams_ptr->k)
		    + "_dt" + std::to_string(simParams_ptr->dt)
		    + "_" + std::to_string(simData_ptr->cells) + "_" + std::to_string(dimensions.width) 
		    + "_" + std::to_string(dimensions.height) + "_" + std::to_string(simData_ptr->stepsRan) 
		    + "_" + std::to_string(simData_ptr->initialSeed);
}

//TODO: an actual reasonable save system : p
//NOTE: values are clamped to [0, 1] for the .pgm image
bool PFM::SimulationControl::saveFieldData(bool savePGM, bool saveBIN, bool saveDAT) const {

	if(!savePGM && !saveBIN && !saveDAT) { return true; } //sucessfully did nothing : )

	if (!m_hasInitialized || m_activeBaseField_ptr == NULL) { return false; }
	
	if(!m_activeBaseField_ptr->isInitialized() || !m_activeBaseField_ptr->hasAllocated()) {
		return false;
	}

	auto dimensions = m_activeBaseField_ptr->getFieldDimensions();

	std::string baseFilename = getFileName(m_lastSimData.stepsRan, false);
	
	int maxColor = 255;
	FILE* fp_pgm = nullptr;
	if(savePGM) {
		fp_pgm = fopen((baseFilename + ".pgm").c_str(), "wb");
		if(fp_pgm == nullptr) { return false; }
		fprintf(fp_pgm, "P5\n%d %d\n%d\n", (int)dimensions.width, (int)dimensions.height, maxColor);
	}

	FILE* fp_bin = nullptr;
	if(saveBIN) {
		fp_bin = fopen((baseFilename + ".bin").c_str(), "wb");
		if(fp_bin == nullptr) { return false; }
	}

	double value;
	for (int j = 0; j < (int)dimensions.height; j++) {
		for (int i = 0; i < (int)dimensions.width; i++) {
		  value = m_activeBaseField_ptr->getDataPoint({i,j});
		  if(saveBIN) { fprintf(fp_bin, "%f ", value); }
		  if(savePGM) { fprintf(fp_pgm, "%c", (char)(std::clamp(value, 0.0, 1.0) * maxColor)); }
		}
	}
	if(saveBIN) { fclose(fp_bin); }
	if(savePGM) { fclose(fp_pgm); }

	FILE* fp_dat = nullptr;
	if(saveDAT) {
		fp_dat = fopen((baseFilename + ".dat").c_str(), "w");
		if(fp_dat == NULL) { return false; }

		fprintf(fp_dat, "%s\n", getSimDataString().c_str());
		fprintf(fp_dat, "%s\n\n", getSimParamsString().c_str());

		size_t numberChecks = getActiveFieldsCheckVectorElements();
		for (size_t i = 0; i < numberChecks; i++) {
			fprintf(fp_dat, "%s\n\n", getActiveFieldsChecksString(i).c_str());
		}
		
		fclose(fp_dat);
	}
	
	return true;
}

bool PFM::SimulationControl::saveFieldData() const {
	return saveFieldData(m_savePGMonIntermediateChecks, m_saveBINonIntermediateChecks, 
		                                                m_saveDATonIntermediateChecks);
}

void PFM::SimulationControl::updateGammaLambda() {
	if(m_lastSimParameters.A == 0) { return; }
	m_lastSimParameters.gamma = std::sqrt(m_lastSimParameters.A * m_lastSimParameters.k) / 6.0;
	m_lastSimParameters.lambda = 2*std::sqrt(m_lastSimParameters.k / m_lastSimParameters.A);
}

void PFM::SimulationControl::updateKandA() {
	if(m_lastSimParameters.lambda == 0) { return; }
	m_lastSimParameters.A = 12 * m_lastSimParameters.gamma / m_lastSimParameters.lambda;
	m_lastSimParameters.k = 3 * m_lastSimParameters.gamma * m_lastSimParameters.lambda;
}

void PFM::SimulationControl::setAused(double newA) {
	m_lastSimParameters.A = newA;

	updateGammaLambda();
}

void PFM::SimulationControl::setKused(double newK) {
	m_lastSimParameters.k = newK;

	updateGammaLambda();
}

void PFM::SimulationControl::setLambdaUsed(double newLambda) {
	m_lastSimParameters.lambda = newLambda;

	updateKandA();
}

void PFM::SimulationControl::setGammaUsed(double newGamma) {
	m_lastSimParameters.gamma = newGamma;

	updateKandA();
}

void PFM::SimulationControl::updatePhysicalParametersFromInternals() {
	updateGammaLambda();
}

void PFM::SimulationControl::setDTused(double newDT) {
	m_lastSimParameters.dt = newDT;
}

void PFM::SimulationControl::setMaxStepsPerCheckAdded(size_t newStepsPerCheckSaved) {
	m_stepsPerCheckSaved = newStepsPerCheckSaved;
}

void PFM::SimulationControl::setMaxTotalChangePerElementPerCheckAdded(double newMaxTotalChangePerCheck) {
	m_absoluteChangePerCheckSaved = std::max(0.0, newMaxTotalChangePerCheck);
}

void PFM::SimulationControl::setIntermediateDATsaves(bool shouldSave) {
	m_saveDATonIntermediateChecks = shouldSave;
}

void PFM::SimulationControl::setIntermediatePGMsaves(bool shouldSave) {
	m_savePGMonIntermediateChecks = shouldSave;
}

void PFM::SimulationControl::setIntermediateBINsaves(bool shouldSave)  {
	m_saveBINonIntermediateChecks = shouldSave;
}

bool PFM::SimulationControl::checkIfShouldStop() {
	if(m_shouldStop || (m_stepsToRun > 0 && m_lastSimData.stepsRan >= m_stepsToRun)) {
		m_shouldStop = true;
	}

	return m_shouldStop;
}

int PFM::SimulationControl::getNumberCells() const {
	return m_lastSimData.cells;
}

double PFM::SimulationControl::getLastCellSeedValue() const {
	return m_lastSimData.cellSeedValue;
}

bool PFM::SimulationControl::shouldStillExpandSeeds() const {
	return m_seedsNeedExpanding;
}

int PFM::SimulationControl::getStepsPerCheckSaved() const {
	return m_stepsPerCheckSaved;
}

double PFM::SimulationControl::getAbsoluteChangePerCheckSaved() const {
	return m_absoluteChangePerCheckSaved;
}

void PFM::SimulationControl::runForSteps(int steps, double lambda, double gamma, double dt,
	                                     PFM::simFuncEnum simulationToRun, 
	                                     PFM::integrationMethods method) {
	if(controller.isSimulationRunning()) { return; }
	if((int)simulationToRun >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS) { return; }
	if((int)method >=  (int)PFM::integrationMethods::TOTAL_METHODS) { return; }

	m_stepsToRun = steps;
	m_lastSimData.simulFunc = simulationToRun;
	m_lastSimData.method = method;
	m_lastSimParameters.lambda = lambda;
	m_lastSimParameters.gamma = gamma;
	updateKandA();
	m_lastSimParameters.dt = dt;
	m_isRunning = true;

	m_stepsThread = std::thread(*(simFunctionsPtrs_arr[(int)simulationToRun]), this, &m_lastSimData.stepsRan, 
		                                                                                &m_isRunning, method);
	m_stepsThread.detach();
}

const simConfig_t* PFM::SimulationControl::getLastSimConfigPtr() const {
	return &m_lastSimData;
}

std::string PFM::SimulationControl::getSimDataString() const { 
	return m_lastSimData.getSimDataString();
}

simParameters_t* PFM::SimulationControl::getLastSimParametersPtr() {
	return &m_lastSimParameters;
}

std::string PFM::SimulationControl::getSimParamsString() const { 
	return m_lastSimParameters.getSimParamsString();
}

size_t PFM::SimulationControl::getActiveFieldsCheckVectorElements() const {
	return getActiveFieldConstPtr()->getCheckVectorConstPtr()->size();
}

std::string PFM::SimulationControl::getActiveFieldsChecksString(size_t checkNumber) const {
	auto checksVec_ptr = getActiveFieldConstPtr()->getCheckVectorConstPtr();
	size_t elements = checksVec_ptr->size();
	if(checkNumber >= elements) { return std::string(""); }
	
	return checksVec_ptr->at(checkNumber).getChecksStr();
}

void PFM::SimulationControl::printSimDataAndParams() const {
	printf("%s\n%s\n", getSimDataString().c_str(), getSimParamsString().c_str());
}

bool PFM::SimulationControl::isSimulationRunning() const {
	return m_isRunning;
}

int PFM::SimulationControl::stepsAlreadyRan() const {
	return m_lastSimData.stepsRan;
}

void PFM::SimulationControl::resetStepsAlreadyRan() {
	m_lastSimData.stepsRan = 0;
}

///These API calls are really just wrappers to calls to methods of the controller:

const PeriodicDoublesLattice2D* PFM::initializeSimulation(simConfig_t config) {
	
	PFM::fieldDimensions_t dimensions;
	dimensions.width = config.width;
	dimensions.height = config.height;

	controller.reinitializeController(dimensions, config.cells, config.initialContidion, 
		                              config.perCellLayer, config.bias, config.initialSeed);
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

bool PFM::saveFieldData(bool savePGM, bool saveBIN, bool saveDAT) {
	return controller.saveFieldData(savePGM, saveBIN, saveDAT);
}
	
void PFM::runForSteps(int stepsToRun, simParameters_t parameters, simConfig_t config) {

	controller.runForSteps(stepsToRun, parameters.lambda, parameters.gamma, parameters.dt, 
		                                      config.simulFunc, config.method);
}

PeriodicDoublesLattice2D* PFM::getActiveFieldPtr() {
	return controller.getActiveFieldPtr();
}

const PeriodicDoublesLattice2D* PFM::getActiveFieldConstPtr() {
	return (const PeriodicDoublesLattice2D*)controller.getActiveFieldPtr();
}

PFM::checkData_t* PFM::getCheckDataPtr() {
	if (!controller.isInitialized() && controller.getActiveFieldPtr()->isInitialized()) { return nullptr; }
	return &controller.getActiveFieldPtr()->checks;
}

const PFM::simConfig_t* PFM::getSimConfigPtr() {
	if (!controller.isInitialized()) { return nullptr; }
	return controller.getLastSimConfigPtr();
}

PFM::simParameters_t* PFM::getSimParamsPtr() {
	if (!controller.isInitialized()) { return nullptr; }
	return controller.getLastSimParametersPtr();
}

void PFM::updatePhysicalParameters() {
	if (!controller.isInitialized()) { return; }
	controller.updatePhysicalParametersFromInternals();
}

void PFM::setMaxStepsPerCheckAdded(size_t newMaxStepsPerCheck) {

}

void PFM::setMaxTotalChangePerElementPerCheckAdded(double newMaxTotalChangePerCheck) {

}

void PFM::setIntermediateDATsaves(bool shouldSave) {

}

void PFM::setIntermediatePGMsaves(bool shouldSave) {

}

void PFM::setIntermediateBINsaves(bool shouldSave) {

}