#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/timeHelpers.hpp"
#include "fAux/API/miscDefines.hpp"
#include "fAux/API/prng.hpp"

#include "PFM_API.hpp"
#include "PFM_tests.hpp"
#include "PFM_defaults.hpp"
#include "simulControl.hpp"

#include <filesystem>

using namespace PFM;

static SimulationControl g_controller;
static uint32_t g_callerKey = PFM::defaultCallerKey;

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
void PFM::SimulationControl::reinitializeController(PFM::simConfig_t config) {
	
	if(isSimulationRunning()) { stop(); }
	releaseFields();

	m_simConfigs = config;
	PFM::fieldDimensions_t dimensions = { (size_t)m_simConfigs.width, (size_t)m_simConfigs.height};
	size_t elements = dimensions.totalElements();
	
	m_lastDphisAndTempKsField_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
			new PeriodicDoublesLattice2D(dimensions, ALL_CELLS_ID)
	);

	std::vector<double> initialData;

	initialData.reserve(elements);

	//The data to be used on initialization:
	int cellSpacing = elements;
	if(m_simConfigs.cells > 1) { cellSpacing /= m_simConfigs.cells; }
	bool centerSingleCell = false; //TODO: pass as parameter

	if(!m_simConfigs.perCellLayer) {
		if (m_simConfigs.initialContidion == PFM::initialConditions::EVENLY_SPACED_INDEX) {
			for (size_t i = 0; i < elements; i++) {
				//To "seed" numberCells equally spaced cells with initialValue, and all others with zero:
				if (m_simConfigs.cells == 1 && centerSingleCell) {
					initialData.push_back(m_simConfigs.cellSeedValue * ( i == (elements + dimensions.width)/2) );
				}
				else {
					initialData.push_back(m_simConfigs.cellSeedValue * ( (i/cellSpacing) == (i/(double)cellSpacing) ) );
				}
			}
			m_seedsNeedExpanding = true;
		}
		else {
			uint64_t seed = m_simConfigs.initialSeed;
			for (size_t i = 0; i < elements; i++) {
				//Each elemnt starts with a random value between -0.5 and 1.5:
				double value = m_simConfigs.bias -0.5 + 2*(AZ::draw1spcg32(&seed)/(double)UINT32_MAX);
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

		m_perCellLaticePtrs.reserve(m_simConfigs.cells);
		for (uint32_t cell = 0; cell < m_simConfigs.cells; cell++) {
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

	m_shouldStop = false;
	setBaseAsActive();
	m_simConfigs.stepsRan = 0;

	m_hasInitialized = true;
}

void PFM::SimulationControl::nonBlockingStop() {
	m_shouldStop = true;
	m_shouldBePaused = false;
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

	return m_simConfigs.stepsRan;
}

void PFM::SimulationControl::pause() {
	m_shouldBePaused = true;
}

void PFM::SimulationControl::resume() {
	m_shouldBePaused = false;
}

bool PFM::SimulationControl::shouldBePaused() const {
	return m_shouldBePaused;
}

void PFM::SimulationControl::mirrorRotatingOnBase() {
	m_baseLattice_ptr.get()->mirrorAllDataFrom(m_rotatingBaseLattice_ptr.get()->getPointerToCurrent());
}

void PFM::SimulationControl::mirrorBaseOnRotating() {
	m_rotatingBaseLattice_ptr.get()->getPointerToCurrent()->mirrorAllDataFrom(m_baseLattice_ptr.get());
}

std::string PFM::getDirAndFileName(int steps, bool calledFromGUI) {
	
	//gather the relevant data:
	auto dimensions = g_controller.getActiveFieldPtr()->getFieldDimensions();
	const PFM::simConfig_t* simConfig_ptr = g_controller.getLastSimConfigPtr();
	const PFM::simParameters_t* simParams_ptr = g_controller.getLastSimParametersPtr();

	std::string dirname = "Sim" + std::to_string((int)simConfig_ptr->simulFunc) 
						+ "_w" + std::to_string(dimensions.width) 
						+ "_h" + std::to_string(dimensions.height)
						+ "_c" + std::to_string(simConfig_ptr->cells)
						+ "_ini" + std::to_string((int)simConfig_ptr->initialContidion) 		    
						+ "_b" + std::to_string(simConfig_ptr->bias)
						+ "m" + std::to_string((int)simConfig_ptr->method)
						+ "_t" + std::to_string(simConfig_ptr->reducedSecondsSinceEpochOnSimCall());	
	
	bool createdDirectory = std::filesystem::create_directory(dirname);
	uint32_t callerKey = PFM::getCallerKey();

	if (createdDirectory && (callerKey != PFM::defaultCallerKey)) {
		//We should let the caller know the path to the new folder:
		FILE* caller_fp = fopen( (std::to_string(callerKey) + ".txt").c_str(), "w");
		if (caller_fp == NULL) {
			//Poor caller
			assert(false && "Couldn't generate the file to let the caller know the path to the new data folder");
		}
		else {
			fprintf(caller_fp, "%s", dirname.c_str());
			fclose(caller_fp);
		}
	}

	std::string baseFilename = "";

	//TODO: do I actually want to mark manual saves?
	//if(calledFromGUI) { fileName += "m"; }

	baseFilename += "s" + std::to_string(simConfig_ptr->stepsRan)
		         + "_A" + std::to_string(simParams_ptr->A) + "_k" + std::to_string(simParams_ptr->k)
		         + "_dt" + std::to_string(simParams_ptr->dt);
		    
	return dirname + "\\" + baseFilename;
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

	std::string baseFilename = getDirAndFileName(m_simConfigs.stepsRan, false);
	
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
		  if(saveBIN) { fwrite(&value, sizeof(value), 1, fp_bin); }
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
			if(m_saveOnDATtheParamsBeforeEachCheck) { 
				fprintf(fp_dat, "%s", getActiveFieldsParamStringBeforeAGivenCheck(i).c_str());
			}
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
	if(m_simParameters.A == 0) { return; }
	m_simParameters.gamma = PFM::getGammaFromKandA(m_simParameters.k, m_simParameters.A);
	m_simParameters.lambda = PFM::getLambdaFromKandA(m_simParameters.k, m_simParameters.A);
}

void PFM::SimulationControl::updateKandA() {
	if(m_simParameters.lambda == 0) { return; }
	m_simParameters.A = 12 * m_simParameters.gamma / m_simParameters.lambda;
	m_simParameters.k = 3 * m_simParameters.gamma * m_simParameters.lambda;
}

void PFM::SimulationControl::setAused(double newA) {
	m_simParameters.A = newA;

	updateGammaLambda();
}

void PFM::SimulationControl::setKused(double newK) {
	m_simParameters.k = newK;

	updateGammaLambda();
}

void PFM::SimulationControl::setLambdaUsed(double newLambda) {
	m_simParameters.lambda = newLambda;

	updateKandA();
}

void PFM::SimulationControl::setGammaUsed(double newGamma) {
	m_simParameters.gamma = newGamma;

	updateKandA();
}

void PFM::SimulationControl::updatePhysicalParametersFromInternals() {
	updateGammaLambda();
}

void PFM::SimulationControl::setDTused(double newDt) {
	m_simParameters.dt = newDt;
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

void PFM::SimulationControl::setSavingOnDATofTheParamsBeforeEachCheck(bool shouldSave)  {
	m_saveOnDATtheParamsBeforeEachCheck = shouldSave;
}

bool PFM::SimulationControl::checkIfShouldStop() {
	if(m_shouldStop || (m_stepsToRun > 0 && m_simConfigs.stepsRan >= m_stepsToRun)) {
		m_shouldStop = true;
	}

	return m_shouldStop;
}

int PFM::SimulationControl::getNumberCells() const {
	return m_simConfigs.cells;
}

double PFM::SimulationControl::getLastCellSeedValue() const {
	return m_simConfigs.cellSeedValue;
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

void PFM::SimulationControl::runForSteps(uint64_t steps, simParameters_t parameters,
	                                     PFM::simFuncEnum simulationToRun, PFM::integrationMethods method) {

	if(g_controller.isSimulationRunning()) { return; }
	if((int)simulationToRun >= (int)PFM::simFuncEnum::TOTAL_SIM_FUNCS) { return; }
	if((int)method >=  (int)PFM::integrationMethods::TOTAL_METHODS) { return; }

	m_stepsToRun = steps;
	m_simConfigs.simulFunc = simulationToRun;
	m_simConfigs.method = method;
	m_simParameters.lambda = parameters.lambda;
	m_simParameters.gamma = parameters.gamma;
	updateKandA();
	m_simParameters.dt = parameters.dt;
	m_simParameters.adaptativeDt = parameters.adaptativeDt;
	m_isRunning = true;

	m_shouldBePaused = m_simConfigs.startPaused;

	m_stepsThread = std::thread(*(simFunctionsPtrs_arr[(int)simulationToRun]), this, &m_simConfigs.stepsRan, 
	                                                      g_controller.getIsPaused_ptr(), &m_isRunning, method);
	m_stepsThread.detach();
}

const simConfig_t* PFM::SimulationControl::getLastSimConfigPtr() const {
	return &m_simConfigs;
}

std::string PFM::SimulationControl::getSimDataString() const { 
	return m_simConfigs.getSimDataString();
}

simParameters_t* PFM::SimulationControl::getLastSimParametersPtr() {
	return &m_simParameters;
}

std::string PFM::SimulationControl::getSimParamsString() const { 
	return m_simParameters.getSimParamsString();
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

std::string PFM::SimulationControl::getActiveFieldsParamStringBeforeAGivenCheck(size_t checkNumber) const {
	auto checksVec_ptr = getActiveFieldConstPtr()->getCheckVectorConstPtr();
	size_t elements = checksVec_ptr->size();
	if(checkNumber >= elements) { return std::string(""); }
	
	return checksVec_ptr->at(checkNumber).getParametersLastUpdateStr();
}

void PFM::SimulationControl::printSimDataAndParams() const {
	printf("%s\n%s\n", getSimDataString().c_str(), getSimParamsString().c_str());
}

bool PFM::SimulationControl::isSimulationRunning() const {
	return m_isRunning;
}

const bool* PFM::SimulationControl::getIsPaused_ptr() const {
	return &m_shouldBePaused;
}

uint64_t PFM::SimulationControl::stepsAlreadyRan() const {
	return m_simConfigs.stepsRan;
}

void PFM::SimulationControl::resetStepsAlreadyRan() {
	m_simConfigs.stepsRan = 0;
}

void PFM::SimulationControl::updateEpochTimeSimCall() {
	m_simConfigs.epochTimeSimCall = std::chrono::system_clock::now();
}

///These API calls are really just wrappers to calls to methods of the controller:

const PeriodicDoublesLattice2D* PFM::initializeSimulation(simConfig_t config) {
	if(isSimulationRunning()) { return g_controller.getActiveFieldPtr(); }

	g_controller.reinitializeController(config);
	return g_controller.getActiveFieldPtr();
}

bool PFM::isSimulationRunning() {
	return g_controller.isSimulationRunning();
}

bool PFM::isSimulationPaused() {
	return PFM::isSimulationRunning() && g_controller.shouldBePaused();
}

int PFM::stopSimulation() {
	return g_controller.stop();
}

void PFM::pauseSimulation() {
	g_controller.pause();
}

void PFM::resumeSimulation() {
	g_controller.resume();
}

uint64_t PFM::getStepsRan() {
	return g_controller.stepsAlreadyRan();
}

void PFM::resetStepsRan() {
	g_controller.resetStepsAlreadyRan();
}

bool PFM::saveFieldData(bool savePGM, bool saveBIN, bool saveDAT) {
	return g_controller.saveFieldData(savePGM, saveBIN, saveDAT);
}

bool PFM::saveFieldDataAccordingToController() {
	return g_controller.saveFieldData();
}

	
void PFM::runForSteps(uint64_t stepsToRun, simParameters_t parameters, simConfig_t config) {

	g_controller.updateEpochTimeSimCall();
	g_controller.runForSteps(stepsToRun, parameters, config.simulFunc, config.method);
}

PeriodicDoublesLattice2D* PFM::getActiveFieldPtr() {
	return g_controller.getActiveFieldPtr();
}

const PeriodicDoublesLattice2D* PFM::getActiveFieldConstPtr() {
	return (const PeriodicDoublesLattice2D*)g_controller.getActiveFieldPtr();
}

PFM::checkData_t* PFM::getCheckDataPtr() {
	if (!g_controller.isInitialized() && g_controller.getActiveFieldPtr()->isInitialized()) { return nullptr; }
	return &g_controller.getActiveFieldPtr()->checks;
}

const PFM::simConfig_t* PFM::getSimConfigPtr() {
	if (!g_controller.isInitialized()) { return nullptr; }
	return g_controller.getLastSimConfigPtr();
}

PFM::simParameters_t* PFM::getSimParamsPtr() {
	if (!g_controller.isInitialized()) { return nullptr; }
	return g_controller.getLastSimParametersPtr();
}

void PFM::updatePhysicalParameters() {
	if (!g_controller.isInitialized()) { return; }
	g_controller.updatePhysicalParametersFromInternals();
}

double PFM::getLambdaFromKandA(double k, double A) {
	return 2*std::sqrt(k / A);
}

double PFM::getGammaFromKandA(double k, double A) {
	return std::sqrt(A * k) / 6.0;
}

PFM::parameterBounds_t PFM::calculateParameterBounds(double k, double A, double dt, uint64_t steps) {
	parameterBounds_t bounds;

	double completelyArbitraryMaxKplusATimesDt = 1.1;
	double completelyArbitraryMaxDtRatioAtStepZero = 0.1;

	//These will be used to find the maximum allowed value of each parameter
	double minValuesForLimitCalc = 0.000001;
	double dividendForLimitCalc = std::max(minValuesForLimitCalc, dt);

	bounds.maxK = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc - A;

	bounds.maxA = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc - k;

	dividendForLimitCalc = std::max(minValuesForLimitCalc, A + k);
	bounds.maxDt = completelyArbitraryMaxKplusATimesDt/dividendForLimitCalc;

	uint64_t stepsRemainingToUnlockDt = PFM::completelyArbitraryStepToUnlockFullDt - steps;
	if(steps > PFM::completelyArbitraryStepToUnlockFullDt) { stepsRemainingToUnlockDt = 0; }
	double fractionOfArbitraryStepsRemaining = 
						  (double)stepsRemainingToUnlockDt/PFM::completelyArbitraryStepToUnlockFullDt;

	bounds.maxDt *= 1.0 - (fractionOfArbitraryStepsRemaining * (1 - completelyArbitraryMaxDtRatioAtStepZero));
	assert(bounds.maxDt > 0 && "Max dt should always be larger than zero");
	
	return bounds;
}

double PFM::calculateMaxAdaptativeDt(const simParameters_t* parameters_ptr, const simConfig_t* config_ptr, const checkData_t* lastCheck_ptr) {
	
	return PFM::calculateParameterBounds(parameters_ptr->k, parameters_ptr->A, parameters_ptr->dt, lastCheck_ptr->step).maxDt;
}

void PFM::setMaxStepsPerCheckAdded(size_t newMaxStepsPerCheck) {
	g_controller.setMaxStepsPerCheckAdded(newMaxStepsPerCheck);
}

void PFM::setMaxTotalChangePerElementPerCheckAdded(double newMaxTotalChangePerCheck) {
	g_controller.setMaxTotalChangePerElementPerCheckAdded(newMaxTotalChangePerCheck);
}

size_t PFM::getMaxStepsPerCheckAdded() {
	return g_controller.getStepsPerCheckSaved();
}

double PFM::getMaxTotalChangePerElementPerCheckAdded() {
	return g_controller.getAbsoluteChangePerCheckSaved();
}

void PFM::setIntermediateDATsaves(bool shouldSave) {
	g_controller.setIntermediateDATsaves(shouldSave);
}

void PFM::setIntermediatePGMsaves(bool shouldSave) {
	g_controller.setIntermediatePGMsaves(shouldSave);
}

void PFM::setIntermediateBINsaves(bool shouldSave) {
	g_controller.setIntermediateBINsaves(shouldSave);
}

void PFM::setSavingOnDATofTheParamsBeforeEachCheck(bool shouldSave) {
	g_controller.setSavingOnDATofTheParamsBeforeEachCheck(shouldSave);
}

void PFM::setCallerKey(uint32_t key) {
	g_callerKey = key;
}

uint32_t PFM::getCallerKey() {
	return g_callerKey;
}