#include "simulControl.hpp"
#include "derivatives.hpp"
#include <assert.h>

//TODO: general reorganization of this mess of a file : )

//TODO: version which accepts the controller and deals with the data (eg, currentAndLast, multi-layers, etc)
void expandCells(PFM::PeriodicDoublesLattice2D* lattice_ptr, float cellRadius, 
	                                 double cellSeedValue, bool invertValueOn);

//"f", the time derivative for numerical integration, for Cahn-Hiliard
inline double chNumericalF(const PFM::neighborhood9_t* neigh_ptr, double k, double A) {
	double phi = neigh_ptr->getCenter();
	double laplacian = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

	return k*laplacian - A*phi*(1-phi)*(1-2*phi);
}

//TODO: implement methods on simulContro.hpp TODO and then these integrations

void explicitEulerCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                          double dt, double chK, double chA, PFM::checkData_t* checks_ptr) {
	
	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	int height = currentStepField_ptr->getFieldDimensions().height;
	int width = currentStepField_ptr->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double dPhi;
	double phi;

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = lastStepField_ptr->getNeighborhood(centerPoint);
		
			dPhi = chNumericalF(&neigh, chK, chA);
			phi = neigh.getCenter() + dt * dPhi;

			currentStepField_ptr->writeDataPoint(centerPoint, phi);

			checks_ptr->density += phi;
			checks_ptr->absoluteChange += std::abs(dPhi);
		}
	}
}

void implicitEulerCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                                            PFM::PeriodicDoublesLattice2D* baseField_ptr, double dt,
	                                               double chK, double chA, PFM::checkData_t* checks_ptr) {

	if(steps < 1) { return; }
	if(steps == 1) { explicitEulerCahnHiliard(rotatingField_ptr, dt, chK, chA, checks_ptr); return; }

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	int height = currentStepField_ptr->getFieldDimensions().height;
	int width = currentStepField_ptr->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double dPhi;
	double phi;
	double phi0;

	for(int step = 0; step < steps; step++) {

		for (int j = 0; j < height; j++) {
				centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
		
				phi0 = baseField_ptr->getDataPoint(centerPoint);

				neigh = lastStepField_ptr->getNeighborhood(centerPoint);
		
				dPhi = chNumericalF(&neigh, chK, chA);
				phi = neigh.getCenter() + dt * dPhi;

				currentStepField_ptr->writeDataPoint(centerPoint, phi);
			}
		}
		rotatingField_ptr->rotatePointers();
		currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
		lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	}

	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;

			phi0 = baseField_ptr->getDataPoint(centerPoint);
			phi = lastStepField_ptr->getDataPoint(centerPoint);

			checks_ptr->density += phi;
			checks_ptr->absoluteChange += std::abs((phi-phi0)/dt);
		}
	}
}

void heunCahnHiliard(PFM::PeriodicDoublesLattice2D* field_ptr,
	                 PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr) {

}

double rungeKutaKnFinalCoef(int steps, int n) {
	return 0;
}

double rungeKutaKnIntermediateCoef(int steps, int n) {
	return 0;
}

void rungeKutaCahnHiliardFirstStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               double chK, double chA, double dt) {

}

void rungeKutaCahnHiliardIntermediateStep(int RKstep, double coefKnFinal, double coefKnInterm,
	                                      int height, int width, 
	                                      PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
										  PFM::PeriodicDoublesLattice2D* field_ptr,
	                                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                                      double chK, double chA, double dt) {

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double phi0, kn;

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				if(RKstep == 0) tempKs_ptr->writeDataPoint(centerPoint, 0); //clear tempKs

				neigh = lastStepField_ptr->getNeighborhood(centerPoint);

				phi0 = neigh.getCenter();
				kn = dt * chNumericalF(&neigh, chK, chA);

				tempKs_ptr->incrementDataPoint(centerPoint, kn * coefKnFinal);

				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + coefKnInterm * kn);
			}
	}
	rotatingField_ptr->rotatePointers();
}

void rungeKutaCahnHiliardFinalStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	                               PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               double chK, double chA, double dt, PFM::checkData_t* checks_ptr) {

}

void rungeKutaCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
						  PFM::PeriodicDoublesLattice2D* field_ptr,
	                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                      double chK, double chA, double dt, PFM::checkData_t* checks_ptr) {

	if(steps < 1) { return; }
	if(steps == 1) { explicitEulerCahnHiliard(rotatingField_ptr, dt, chK, chA, checks_ptr); return;}

	int height = field_ptr->getFieldDimensions().height;
	int width = field_ptr->getFieldDimensions().width;

	double coefKnFinal = rungeKutaKnFinalCoef(steps, 1);
	double coefKnInterm = rungeKutaKnIntermediateCoef(steps, 1);

	rungeKutaCahnHiliardFirstStep(coefKnFinal, coefKnInterm, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, chK, chA, dt);
	if (steps > 2) {
		for(int i = 0; i < (steps - 2); i++) {
			coefKnFinal = rungeKutaKnFinalCoef(steps, i+2);
			coefKnInterm = rungeKutaKnIntermediateCoef(steps, i+2);

			rungeKutaCahnHiliardIntermediateStep(i, coefKnFinal, coefKnInterm, height, width, 
				                                 rotatingField_ptr, field_ptr, tempKs_ptr, chK, chA, dt);
		}
	}

	coefKnFinal = rungeKutaKnFinalCoef(steps, steps);
	coefKnInterm = rungeKutaKnIntermediateCoef(steps, steps);

	rungeKutaCahnHiliardFinalStep(coefKnFinal, coefKnInterm, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, chK, chA, dt, checks_ptr);
}

//Runs Chan-Hiliard on one layer per cell and adds the results into a base layer
void PFM::multiLayerCHsim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//TODO: Implement : )
	puts("multi-layer simulation not it implemented");
	*isRunning_ptr = false;
}

//TODO: separate initialization, simulation and wrap-up on the following functions

//Runs Chan-Hiliard on a single base layer. Keeps current and last step and uses both to calculate change
void PFM::singleLayerCHsimCurrentAndOld_fn(SimulationControl* controller_ptr, int* stepCount_ptr, 
	                                                                         bool* isRunning_ptr) {

	//TODO: actually rotate the current and old data layers : p

	const bool invertField = true;
	const double k = 1;
	const double A = 0.5;
	double expectedInterfaceWidth = std::sqrt(2*k/A);
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1/std::sqrt(2);
	const int stepsPerCheckSaved = controller_ptr->getStepsPerCheckSaved();

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto rotBaseField_ptr = controller_ptr->getRotatingBaseFieldPtr();
	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
	auto currentStepField_ptr = rotBaseField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotBaseField_ptr->getPointerToLast();
	
	//auto layerFieldsVectorPtr = controller_ptr->getLayerFieldsVectorPtr();

	const int width = (int)baseField_ptr->getFieldDimensions().width;
	const int height = (int)baseField_ptr->getFieldDimensions().height;

	//CELL EXPANSION (EXTRACT):

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

	//ACTUAL STEPS:

	PFM::coordinate_t centerPoint;
	double phi;
	double phi0;
	double kn;
	double dPhi;
	neighborhood9_t neigh;
	PFM::checkData_t checkData;

	auto tempKsField_ptr = controller_ptr->getLastDphiFieldPtr();

	while(!controller_ptr->checkIfShouldStop()) {

		checkData.density = 0;
		checkData.absoluteChange = 0;
		checkData.step = *stepCount_ptr;

		//RK4:
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				tempKsField_ptr->writeDataPoint(centerPoint, 0); //clear tempKs

				neigh = lastStepField_ptr->getNeighborhood(centerPoint);
				
				kn = dt * chNumericalF(&neigh, k, A);
				tempKsField_ptr->incrementDataPoint(centerPoint, (1.0/6) * kn);
				
				phi0 = baseField_ptr->getDataPoint(centerPoint);
				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + kn);
			}
		}
		rotBaseField_ptr->rotatePointers();

		currentStepField_ptr = rotBaseField_ptr->getPointerToCurrent();
		lastStepField_ptr = rotBaseField_ptr->getPointerToLast();
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				neigh = lastStepField_ptr->getNeighborhood(centerPoint);
				
				kn = dt * chNumericalF(&neigh, k, A);
				tempKsField_ptr->incrementDataPoint(centerPoint, (1.0/3) * kn);
				
				phi0 = baseField_ptr->getDataPoint(centerPoint);
				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + 0.5 * kn);
			}
		}
		rotBaseField_ptr->rotatePointers();

		currentStepField_ptr = rotBaseField_ptr->getPointerToCurrent();
		lastStepField_ptr = rotBaseField_ptr->getPointerToLast();
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				neigh = lastStepField_ptr->getNeighborhood(centerPoint);
				
				kn = dt * chNumericalF(&neigh, k, A);
				tempKsField_ptr->incrementDataPoint(centerPoint, (1.0/3) * kn);
				
				phi0 = baseField_ptr->getDataPoint(centerPoint);
				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + 0.5 * kn);
			}
		}
		rotBaseField_ptr->rotatePointers();

		currentStepField_ptr = rotBaseField_ptr->getPointerToCurrent();
		lastStepField_ptr = rotBaseField_ptr->getPointerToLast();
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				neigh = lastStepField_ptr->getNeighborhood(centerPoint);

				kn = dt * chNumericalF(&neigh, k, A);
				tempKsField_ptr->incrementDataPoint(centerPoint, (1.0/6) * kn);
				
				phi0 = baseField_ptr->getDataPoint(centerPoint);
				dPhi = tempKsField_ptr->getDataPoint(centerPoint);
				phi = phi0 + dPhi;

				baseField_ptr->writeDataPoint(centerPoint, phi);		

				checkData.density += phi;
				checkData.absoluteChange += std::abs(dPhi/dt);
			}
		}

		int checksPerPrintout = 5;
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			checksPerPrintout = 1;
		#endif
		
		if(*stepCount_ptr % stepsPerCheckSaved == 0) { 
			checkData.density /= (width * height);
			checkData.absoluteChange /= (width * height);
			currentStepField_ptr->addFieldCheckData(checkData);

			if(*stepCount_ptr % (stepsPerCheckSaved * checksPerPrintout) == 0) {
				printf("steps: %d - density: %f - absolute change (last step): %f\n", 
							*stepCount_ptr, checkData.density, checkData.absoluteChange); 
			}
		}	

		*stepCount_ptr += 1;
		rotBaseField_ptr->rotatePointers();
	}

	*isRunning_ptr = false;
}

//Runs Chan-Hiliard on a single base layer. OVERWRITES values mid-step (doesn't keep last step data separatedly)
void PFM::singleLayerCHsimOnlyCurrent_fn(SimulationControl* controller_ptr, int* stepCount_ptr, 
	                                                                       bool* isRunning_ptr) {
	
	const bool invertField = false;
	const uint32_t backwardEulerExtraSteps = 0;
	const double k = 0.5;
	const double A = 1;
	double expectedInterfaceWidth = std::sqrt(2*k/A);
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1/std::sqrt(2);
	const int stepsPerCheckSaved = controller_ptr->getStepsPerCheckSaved();

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto baseField_ptr = controller_ptr->getRotatingBaseFieldPtr();
	auto currentStepField_ptr = baseField_ptr->getPointerToCurrent();

	//auto layerFieldsVectorPtr = controller_ptr->getLayerFieldsVectorPtr();

	const int width = (int)currentStepField_ptr->getFieldDimensions().width;
	const int height = (int)currentStepField_ptr->getFieldDimensions().height;

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt((float)controller_ptr->getNumberCells()));

	if(controller_ptr->shouldStillExpandSeeds()) {
		expandCells(currentStepField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
	}
	else if (invertField) {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				//TODO: make sure this inversion actually works in general
				currentStepField_ptr->writeDataPoint({i,j}, 1 - currentStepField_ptr->getDataPoint({i,j}));
			}
		}
	}
	
	double radiusToWidth = initialCellRadius / expectedInterfaceWidth;
	printf("Radius to expected interface width = %f (radius: %f, width: %f)\n", 
		             radiusToWidth, initialCellRadius, expectedInterfaceWidth );

	PFM::coordinate_t centerPoint;
	double phi;
	double phi0;
	double laplacian;	
	double change;
	double newValue;
	neighborhood9_t neigh;
	PFM::checkData_t checkData;

	while(!controller_ptr->checkIfShouldStop()) {

		checkData.density = 0;
		checkData.absoluteChange = 0;
		checkData.step = *stepCount_ptr;

		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;

				neigh = currentStepField_ptr->getNeighborhood(centerPoint);
				phi0 = neigh.getCenter();

				//for some backward euler steps (may be "off" if backwardEulerExtraSteps = 0) :
				//TODO: tests needed
				for (uint32_t k = 0; k < backwardEulerExtraSteps; k++) {

					laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
					phi = neigh.getCenter();

					neigh.setCenter(phi0 + dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi)) );
				}
				
				laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
				phi = neigh.getCenter();

				change = dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi));
				newValue = phi0 + change;
				
				checkData.density += newValue;
				checkData.absoluteChange += std::abs(change);

				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + change);
			}
		}

		int checksPerPrintout = 5;
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			checksPerPrintout = 1;
		#endif
		
		if(*stepCount_ptr % stepsPerCheckSaved == 0) { 
			checkData.density /= (width * height);
			checkData.absoluteChange /= (width * height);
			currentStepField_ptr->addFieldCheckData(checkData);

			if(*stepCount_ptr % (stepsPerCheckSaved * checksPerPrintout) == 0) {
				printf("steps: %d - density: %f - absolute change (last step): %f\n", 
							*stepCount_ptr, checkData.density, checkData.absoluteChange); 
			}
		}	

		*stepCount_ptr += 1;
	}

	*isRunning_ptr = false;
}

//Test simulation: the pixels initialized as non-zero should "diffuse" up and to the right
//(not really diffuse, more like reinforce - eventually everything should be "maximal" and then loop)
void PFM::dataAndControllerTest_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	
	const double diffusionFactor = 0.025;
	const double maxValue = 1;
	controller_ptr->setKused(diffusionFactor);
	controller_ptr->setAused(maxValue);	
	controller_ptr->setDTused(1);	

	auto field_ptr = controller_ptr->getRotatingBaseFieldPtr()->getPointerToCurrent();

	const int width = (int)field_ptr->getFieldDimensions().width;
	const int height = (int)field_ptr->getFieldDimensions().height;
	
	double value, valueAbove, valueToTheRight;

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
	if(invertValueOn) { valueOn = 0.5 + (0.5 - valueOn); } //reflect around 0.5

	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			double value = lattice_ptr->getDataPoint({i,j});
			if(value < 0) {
				lattice_ptr->writeDataPoint({i,j}, valueOn);
			}
			else {
				lattice_ptr->writeDataPoint({i,j}, 1 - valueOn);
			}
		}
	}
}