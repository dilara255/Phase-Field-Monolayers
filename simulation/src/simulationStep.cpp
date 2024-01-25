#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/miscDefines.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "PFM_dataDefaults.hpp"

#include "simulControl.hpp"
#include "rateOfChangeFunctions.hpp"
#include "numericalIntegration.hpp"
#include "derivatives.hpp"

//TODO: version which accepts the controller and deals with the data (eg, currentAndLast, multi-layers, etc)
void expandCells(PFM::PeriodicDoublesLattice2D* lattice_ptr, float cellRadius, 
	                                 double cellSeedValue, bool invertValueOn);

void updatedAndSaveChecks(PFM::checkData_t* checks_ptr, const uint64_t step, const int stepsPerCheckSaved, 
	                      const double absoluteChangePerElementPerCheckSaved, const double dt);

void preProccessFieldsAndUpdateController(PFM::SimulationControl* controller_ptr, 
	                                      const double initialCellDiameterDensity, 
	                                      const double k, const double A, const double dt,
	                                      const double expectedInterfaceWidth, const bool invertField);

//Runs Chan-Hiliard on one layer per cell and adds the results into a base layer
void PFM::multiLayerCHsim_fn(SimulationControl* controller_ptr, uint64_t* stepCount_ptr, 
	                         const bool* shouldPause_ptr, bool* isRunning_ptr, 
	                         integrationMethods method = integrationMethods::FTCS) {
	//TODO: Implement : )
	puts("multi-layer simulation not it implemented");
	*isRunning_ptr = false;
}

//Runs Chan-Hiliard on a single layer. Keeps current and last step and uses both to calculate change
void PFM::singleLayerCHsim_fn(SimulationControl* controller_ptr, uint64_t* stepCount_ptr, 
							  const bool* shouldPause_ptr, bool* isRunning_ptr, 
							  integrationMethods method = integrationMethods::FTCS) {

	//Parameters for the steps:
	//TODO: extract the parameters
	const bool invertField = false;

	//we're using fixed dx = dy = 1, h2 = 2
	//TODO: add h as parameter
	//static_assert(dt <= 1/(4*k), "dt too large for k"); //Max stable dt for FTCS heat diffusion, used as a ballpark
	
	const double initialCellDiameterDensity = 1.0/2;

	double dt = controller_ptr->getLastSimParametersPtr()->dt;
	double A = controller_ptr->getLastSimParametersPtr()->A;
	double k = controller_ptr->getLastSimParametersPtr()->k;
	double lambda = controller_ptr->getLastSimParametersPtr()->lambda;

	//Preparation to start stepping:
	preProccessFieldsAndUpdateController(controller_ptr, initialCellDiameterDensity, k, A, dt, lambda, invertField);
	
	controller_ptr->setBaseAsActive();
	auto checks_ptr = controller_ptr->getActiveFieldsCheckDataPtr();
	auto rotBaseField_ptr = controller_ptr->getRotatingBaseFieldPtr();
	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
	auto tempKsAndDphis_ptr = controller_ptr->getLastDphisAndTempKsFieldPtr();
	
	controller_ptr->printSimDataAndParams();
	checks_ptr->parametersOnLastCheck = *controller_ptr->getLastSimParametersPtr();
	updatedAndSaveChecks(checks_ptr, *stepCount_ptr, controller_ptr->getStepsPerCheckSaved(), 
		                                controller_ptr->getAbsoluteChangePerCheckSaved(), dt);

	//We call a first dt update to sanitize possibly too high initial dt
	//This is done before the checks and their printing so they reflect the original intent
	controller_ptr->updateDt();

	//Start paused?
	while(*shouldPause_ptr) { AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI)); }

	//The actual steps:
	auto params_ptr = controller_ptr->getLastSimParametersPtr();
	while(!controller_ptr->checkIfShouldStop()) {
	
		//Update dt each step so we can deal with the initial few steps, adaptative dt, and client changes
		controller_ptr->updateDt();
		checks_ptr->clearLastStepsChanges();

		//TODO: maybe pull into an "updateLocalParameters" function?
		dt = params_ptr->dt;
		A = params_ptr->A;
		k = params_ptr->k;

		switch (method)
		{
		case PFM::integrationMethods::FTCS:
			controller_ptr->setBaseAsActive();
			N_INT::TD::CH::ftcsStep(baseField_ptr, tempKsAndDphis_ptr, dt, k, A, checks_ptr);
			break;
		case PFM::integrationMethods::FTCS_WITH_SUBS:
			controller_ptr->setBaseAsActive();
			N_INT::TD::CH::ftcsStepWithSubsteps(baseField_ptr, tempKsAndDphis_ptr, dt, k, A, checks_ptr,
				                                rotBaseField_ptr->getPointerToCurrent());
			break;
		case PFM::integrationMethods::HEUN:
			controller_ptr->setRotatingLastAsActive();
			N_INT::TD::CH::heunStep(rotBaseField_ptr, tempKsAndDphis_ptr, dt, k, A, checks_ptr);
			break;
		case PFM::integrationMethods::RK2:
			controller_ptr->setBaseAsActive();
			N_INT::TD::CH::rungeKuttaStep(N_INT::rungeKuttaOrder::TWO, rotBaseField_ptr, baseField_ptr, 
			                                             tempKsAndDphis_ptr, dt, k, A, checks_ptr);
			break;
		case PFM::integrationMethods::RK4:
			controller_ptr->setBaseAsActive();
			N_INT::TD::CH::rungeKuttaStep(N_INT::rungeKuttaOrder::FOUR, rotBaseField_ptr, baseField_ptr, 
			                                             tempKsAndDphis_ptr, dt, k, A, checks_ptr);
			break;
		default:
			puts("Bad Method Code"); assert(false);
			return;
		}

		*stepCount_ptr += 1;

		checks_ptr->parametersOnLastCheck = *controller_ptr->getLastSimParametersPtr();
		updatedAndSaveChecks(checks_ptr, *stepCount_ptr, controller_ptr->getStepsPerCheckSaved(), 
			                                controller_ptr->getAbsoluteChangePerCheckSaved(), dt);
		
		//Pause?
		while(*shouldPause_ptr) { AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI)); }
	}

	//Done:
	*isRunning_ptr = false;
}

//Test simulation: the pixels initialized as non-zero should "diffuse" up and to the right
//(not really diffuse, more like reinforce - eventually everything should be "maximal" and then loop)
void PFM::dataAndControllerTest_fn(SimulationControl* controller_ptr, uint64_t* stepCount_ptr, 
								   const bool* shouldPause_ptr, bool* isRunning_ptr, 
								   integrationMethods method = integrationMethods::FTCS) {
	
	const double diffusionFactor = 0.025;
	const double maxValue = 1;
	controller_ptr->setKused(diffusionFactor);
	controller_ptr->setAused(maxValue);	
	controller_ptr->setIntendedDt(1);	

	auto field_ptr = controller_ptr->getRotatingBaseFieldPtr()->getPointerToCurrent();

	const int width = (int)field_ptr->getFieldDimensions().width;
	const int height = (int)field_ptr->getFieldDimensions().height;
	
	double value, valueAbove, valueToTheRight;

	//Start paused?
	while(*shouldPause_ptr) { AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI)); }

	while(!controller_ptr->checkIfShouldStop()) {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {

				value = field_ptr->getDataPoint({i,j});
				valueAbove = field_ptr->getDataPoint({i,j-1});
				valueToTheRight = field_ptr->getDataPoint({i+1,j});

				field_ptr->writeDataPoint({i,j-1}, valueAbove + (diffusionFactor * value) );
				field_ptr->writeDataPoint({i+1,j}, valueToTheRight + (diffusionFactor * value) );

				if(value >= maxValue) { field_ptr->writeDataPoint({i,j}, (value - maxValue)); }
			}
		}

		*stepCount_ptr += 1;
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			if(*stepCount_ptr % 100 == 0) { printf("steps: %d\n", *stepCount_ptr); }
		#endif

		//Pause?
		while(*shouldPause_ptr) { AZ::hybridBusySleepForMicros(std::chrono::microseconds(MICROS_IN_A_MILLI)); }
	}
	
	*isRunning_ptr = false;
}

void expandCells(PFM::PeriodicDoublesLattice2D* lattice_ptr, float cellRadius, 
	                                 double cellSeedValue, bool invertValueOn) {
	assert(cellRadius >= 0);

	const int width = (int)lattice_ptr->getFieldDimensions().width;
	const int height = (int)lattice_ptr->getFieldDimensions().height;

	PFM::coordinate_t centerPoint;
	int x, y;
	int testRadius = (int)std::ceil(cellRadius);
	double centerValue;
	const double markerValue = -1;
	const double sqrdRadius = cellRadius * cellRadius;
	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;

			centerValue = lattice_ptr->getDataPoint(centerPoint);
			if (centerValue ==  cellSeedValue) {
				
				for (int j_delta = -testRadius; j_delta <= testRadius; j_delta++) {
					for (int i_delta = -testRadius; i_delta <= testRadius; i_delta++) {
						
						double sqrdDistanceCenter = j_delta * j_delta + i_delta * i_delta;

						if (sqrdDistanceCenter <= sqrdRadius) {
							x = i + i_delta;
							y = j + j_delta;

							lattice_ptr->writeDataPoint({x,y}, markerValue);
						}	
					}
				}
			}
		}
	}

	double valueOn = cellSeedValue;
	double inverseValueOn = 0.5 + (0.5 - valueOn); //reflect around 0.5
	double valueOff = inverseValueOn;
	if(invertValueOn) { valueOff = valueOn; valueOn = inverseValueOn; } 
	double newValue;
	lattice_ptr->checks.zeroOut();
	
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {

			double value = lattice_ptr->getDataPoint({i,j});

			if(value < 0) { newValue = valueOn; }
			else { newValue = valueOff; }
			
			lattice_ptr->writeDataPoint({i,j}, newValue);
			lattice_ptr->checks.lastDensity += newValue;
		}
	}

	lattice_ptr->checks.lastDensity /= (width * height);
}

void preProccessFieldsAndUpdateController(PFM::SimulationControl* controller_ptr, 
	                                      const double initialCellDiameterDensity, 
	                                      const double k, const double A, const double dt,
	                                      const double expectedInterfaceWidth, const bool invertField) {
	
	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setIntendedDt(dt);

	auto rotBaseField_ptr = controller_ptr->getRotatingBaseFieldPtr();
	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
	auto currentStepField_ptr = rotBaseField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotBaseField_ptr->getPointerToLast();
	
	//auto layerFieldsVectorPtr = controller_ptr->getLayerFieldsVectorPtr();

	const int width = (int)baseField_ptr->getFieldDimensions().width;
	const int height = (int)baseField_ptr->getFieldDimensions().height;

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt((float)controller_ptr->getNumberCells()));

	if(controller_ptr->shouldStillExpandSeeds()) {
		expandCells(currentStepField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
		expandCells(lastStepField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
		expandCells(baseField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
	}
	else if (invertField) {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				//TODO: make sure this inversion actually works in general
				currentStepField_ptr->writeDataPoint({i,j}, 1 - currentStepField_ptr->getDataPoint({i,j}));
				lastStepField_ptr->writeDataPoint({i,j}, 1 - lastStepField_ptr->getDataPoint({i,j}));
				baseField_ptr->writeDataPoint({i,j}, 1 - baseField_ptr->getDataPoint({i,j}));
			}
		}
	}
	
	double radiusToWidth = initialCellRadius / expectedInterfaceWidth;
	printf("Radius to expected interface width = %f (radius: %f, width: %f)\n", 
		             radiusToWidth, initialCellRadius, expectedInterfaceWidth );
}

void updatedAndSaveChecks(PFM::checkData_t* checks_ptr, const uint64_t step, const int stepsPerCheckSaved, 
	                      const double absoluteChangePerElementPerCheckSaved, const double dt) {
	
	checks_ptr->step = step;
	checks_ptr->totalTime += dt;

	checks_ptr->absoluteChange += checks_ptr->absoluteChangeLastStep;
	checks_ptr->sumOfsquaresOfChanges += checks_ptr->sumOfsquaresOfChangesLastStep;

	const int elements = PFM::getActiveFieldConstPtr()->getNumberOfActualElements();
	checks_ptr->totalAbsoluteChangeSinceLastSave = 
				checks_ptr->remainingChangeSinceSaveOnLastCheck  + checks_ptr->absoluteChange / elements;
	
	static uint64_t stepsSinceCheckPrintout = 0;
	if(( (step - checks_ptr->stepsAtLastCheck) % stepsPerCheckSaved == 0) || 
		 checks_ptr->totalAbsoluteChangeSinceLastSave >= absoluteChangePerElementPerCheckSaved) { 
		
		checks_ptr->referenceDt = dt;

		//The max is there in case update is called before the simulation begins (to record initial condition)
		checks_ptr->stepsDuringLastCheckPeriod = std::max(1llu, checks_ptr->step - checks_ptr->stepsAtLastCheck);
		checks_ptr->stepsAtLastCheck = checks_ptr->step;
		//The max is there for the same reason
		checks_ptr->timeDuringLastCheckPeriod = std::max(1.0, checks_ptr->totalTime - checks_ptr->timeAtLastcheck);
		checks_ptr->timeAtLastcheck = checks_ptr->totalTime;

		checks_ptr->lastDensityChange = checks_ptr->densityChange;
		checks_ptr->lastDensity += (checks_ptr->lastDensityChange / elements);
		checks_ptr->lastAbsoluteChangePerElement = checks_ptr->absoluteChange / elements;
		
		double changeMeanSquare = checks_ptr->sumOfsquaresOfChanges / elements;
		double squaredMeanChange = (checks_ptr->densityChange / elements) * (checks_ptr->densityChange / elements);
		
		checks_ptr->lastRmsAbsChange = std::sqrt(changeMeanSquare);
		checks_ptr->lastAbsChangeStdDev = std::sqrt(changeMeanSquare - squaredMeanChange);
	
		//How much did the change overshoot the threshold?
		checks_ptr->remainingChangeSinceSaveOnLastCheck = PFM::absoluteChangeRemainingFactor 
				* (checks_ptr->totalAbsoluteChangeSinceLastSave - absoluteChangePerElementPerCheckSaved);

		//In case change was actually *smaller* than the threshold:
		if (checks_ptr->remainingChangeSinceSaveOnLastCheck < 0) {
			//Just keep the change for the next cycle:
			checks_ptr->remainingChangeSinceSaveOnLastCheck = checks_ptr->totalAbsoluteChangeSinceLastSave;
		}

		PFM::getActiveFieldPtr()->addFieldCheckData(*checks_ptr);
		PFM::saveFieldDataAccordingToController();
		
		int checksPerPrintout = 5; //TODO: parametrize this
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			checksPerPrintout = 1;
		#endif
		uint64_t minPrintOutPeriod = stepsPerCheckSaved * checksPerPrintout;
		

		if(stepsSinceCheckPrintout >= minPrintOutPeriod || step == 0) {
			printf("%s", checks_ptr->getChecksStr().c_str());
			stepsSinceCheckPrintout = 0;
		}

		checks_ptr->clearCurrentChanges();
	}	
	stepsSinceCheckPrintout++;
}