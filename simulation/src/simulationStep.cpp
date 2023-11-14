#include "simulControl.hpp"
#include "derivatives.hpp"
#include <assert.h>

//TODO: separate initialization, simulation and wrap-up

//TODO: version which accepts the controller and deals with the data (eg, currentAndLast, multi-layers, etc)
void expandCells(PFM::PeriodicDoublesLattice2D* lattice_ptr, float cellRadius, 
	                                 double cellSeedValue, bool invertValueOn);

//Runs Chan-Hiliard on one layer per cell and adds the results into a base layer
void PFM::multiLayerCHsim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//TODO: Implement : )
	puts("multi-layer simulation not it implemented");
	*isRunning_ptr = false;
}

//Runs Chan-Hiliard on a single base layer. Keeps current and last step and uses both to calculate change
void PFM::singleLayerCHsimCurrentAndOld_fn(SimulationControl* controller_ptr, int* stepCount_ptr, 
	                                                                         bool* isRunning_ptr) {

	//TODO: actually rotate the current and old data layers : p

	const bool invertField = false;
	const uint32_t extraSubsteps = 2;
	const double k = 1;
	const double A = 0.5;
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1/std::sqrt(2);
	const int stepsPerCheckSaved = controller_ptr->getStepsPerCheckSaved();

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
	auto currentStepField_ptr = baseField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = baseField_ptr->getPointerToLast();
	
	//auto layerFieldsVectorPtr = controller_ptr->getLayerFieldsVectorPtr();

	const int width = (int)currentStepField_ptr->getFieldDimensions().width;
	const int height = (int)currentStepField_ptr->getFieldDimensions().height;

	//CELL EXPANSION (EXTRACT):

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt((float)controller_ptr->getNumberCells()));

	if(controller_ptr->shouldStillExpandSeeds()) {
		expandCells(currentStepField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
		expandCells(lastStepField_ptr, initialCellRadius, 
			            controller_ptr->getLastCellSeedValue(), invertField);
	}
	else if (invertField) {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				//TODO: make sure this inversion actually works in general
				currentStepField_ptr->writeDataPoint({i,j}, 1 - currentStepField_ptr->getDataPoint({i,j}));
				lastStepField_ptr->writeDataPoint({i,j}, 1 - lastStepField_ptr->getDataPoint({i,j}));
			}
		}
	}
	
	double expectedInterfaceWidth = std::sqrt(2*k/A);
	double radiusToWidth = initialCellRadius / expectedInterfaceWidth;
	printf("Radius to expected interface width = %f (radius: %f, width: %f)\n", 
		             radiusToWidth, initialCellRadius, expectedInterfaceWidth );

	//ACTUAL STEPS:

	PFM::coordinate_t centerPoint;
	double phi;
	double phi0;
	double laplacian;	
	double change;
	double newValue;
	neighborhood9_t neigh;
	PFM::checkData_t checkData;

	auto lastDphisField_ptr = controller_ptr->getLastDphiFieldPtr();

	while(!controller_ptr->checkIfShouldStop()) {

		checkData.density = 0;
		checkData.absoluteChange = 0;
		checkData.step = *stepCount_ptr;

		currentStepField_ptr = baseField_ptr->getPointerToCurrent();
		lastStepField_ptr = baseField_ptr->getPointerToLast();

		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;

				neigh = lastStepField_ptr->getNeighborhood(centerPoint);
				phi0 = neigh.getCenter();

				for (uint32_t k = 0; k < extraSubsteps; k++) {

					laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
					phi = neigh.getCenter();

					neigh.setCenter(phi0 + dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi)) );
				}
				
				laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
				phi = neigh.getCenter();

				change = 0.5* (dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi))) 
					   + 0.5 * lastDphisField_ptr->getDataPoint(centerPoint);

				lastDphisField_ptr->writeDataPoint(centerPoint, change);

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
		baseField_ptr->rotatePointers();
	}

	*isRunning_ptr = false;
}

//Runs Chan-Hiliard on a single base layer. OVERWRITES values mid-step (doesn't keep last step data separatedly)
void PFM::singleLayerCHsimOnlyCurrent_fn(SimulationControl* controller_ptr, int* stepCount_ptr, 
	                                                                       bool* isRunning_ptr) {
	
		const bool invertField = false;
	const uint32_t extraSubsteps = 2;
	const double k = 1;
	const double A = 0.5;
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1/std::sqrt(2);
	const int stepsPerCheckSaved = controller_ptr->getStepsPerCheckSaved();

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
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
	
	double expectedInterfaceWidth = std::sqrt(2*k/A);
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

				for (uint32_t k = 0; k < extraSubsteps; k++) {

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

	auto field_ptr = controller_ptr->getBaseFieldPtr()->getPointerToCurrent();

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
	if(invertValueOn) valueOn = 1 - valueOn;

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