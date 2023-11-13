#include "simulControl.hpp"
#include "derivatives.hpp"
#include <assert.h>

//TODO: separate initialization, simulation and wrap-up

void expandCells(PFM::PeriodicDoublesLattice2D* lattice_ptr, float cellRadius, 
	                                 double cellSeedValue, bool invertValueOn);

void PFM::multiLayerCHsim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//Runs Chan-Hiliard on one layer per cell and adds the resultos into a base layer

	const bool invertField = false;
	const uint32_t extraSubsteps = 2;
	const double k = 1;
	const double A = 0.5;
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1; //always a bit less because integers. Also, seeds not packed
	//TODO: on the seeding function, define spacing (ie, option to use a spatial configuration for best packing)

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto baseField_ptr = controller_ptr->getBaseFieldPtr();
	auto layerFieldsVectorPtr = controller_ptr->getLayerFieldsVectorPtr();

	const int width = (int)baseField_ptr->getFieldDimensions().width;
	const int height = (int)baseField_ptr->getFieldDimensions().height;

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt(2) * std::sqrt((float)controller_ptr->getNumberCells()));

	expandCells(baseField_ptr, initialCellRadius, controller_ptr->getLastCellSeedValue(), invertField);
	
	double expectedInterfaceWidth = std::sqrt(2*k/A);
	double radiusToWidth = initialCellRadius / expectedInterfaceWidth;
	printf("Radius to expected interface width = %f (radius: %f, width: %f)\n", 
		             radiusToWidth, initialCellRadius, expectedInterfaceWidth );

	PFM::coordinate_t centerPoint;
	double phi;
	double phi0;
	double laplacian;	
	neighborhood9_t neigh;

	while(!controller_ptr->checkIfShouldStop()) {
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;

				neigh = baseField_ptr->getNeighborhood(centerPoint);
				phi0 = neigh.getCenter();

				for (uint32_t k = 0; k < extraSubsteps; k++) {

					laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
					phi = neigh.getCenter();

					neigh.setCenter(phi0 + dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi)) );
				}
				
				laplacian = PFM::laplacian9pointsAroundNeighCenter(&neigh);
				phi = neigh.getCenter();

				baseField_ptr->incrementDataPoint(centerPoint, dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi)) );
			}
		}

		*stepCount_ptr += 1;
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			if(*stepCount_ptr % 1000 == 0) { printf("steps: %d\n", *stepCount_ptr); }
		#else
			if(*stepCount_ptr % 100000 == 0) { printf("steps: %d\n", *stepCount_ptr); }
		#endif
	}

	*isRunning_ptr = false;
}

void PFM::singleLayerCHsim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//Runs Chan-Hiliard on a single layer
	
	const bool invertField = false;
	const double k = 1;
	const double A = 0.125;
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1; //always a bit less because integers. Also, seeds not packed
	//TODO: on the seeding function, define spacing (ie, option to use a spatial configuration for best packing)

	controller_ptr->setKused(k);
	controller_ptr->setAused(A);
	controller_ptr->setDTused(dt);

	auto field_ptr = controller_ptr->getBaseFieldPtr();

	const int width = (int)field_ptr->getFieldDimensions().width;
	const int height = (int)field_ptr->getFieldDimensions().height;

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt(2) * std::sqrt((float)controller_ptr->getNumberCells()));

	expandCells(field_ptr, initialCellRadius, controller_ptr->getLastCellSeedValue(), invertField);
	
	double expectedInterfaceWidth = std::sqrt(2*k/A);
	double radiusToWidth = initialCellRadius / expectedInterfaceWidth;
	printf("Radius to expected interface width = %f (radius: %f, width: %f)\n", 
		             radiusToWidth, initialCellRadius, expectedInterfaceWidth );

	PFM::coordinate_t centerPoint;
	double phi;
	double laplacian;	

	while(!controller_ptr->checkIfShouldStop()) {
		for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;

				phi = field_ptr->getDataPoint(centerPoint);
				laplacian = PFM::laplacian5points(field_ptr, centerPoint);

				field_ptr->incrementDataPoint(centerPoint, dt*(k*laplacian - A*phi*(1-phi)*(1-2*phi)) );
			}
		}

		*stepCount_ptr += 1;
		#ifdef AS_DEBUG //TODO: this is a definition from the build system which should change
			if(*stepCount_ptr % 100 == 0) { printf("steps: %d\n", *stepCount_ptr); }
		#endif
	}
	
	*isRunning_ptr = false;
}

void PFM::dataAndControllerTest_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//Test simulation: the pixels initialized as non-zero should "diffuse" up and to the right
	//(not really diffuse, more like reinforce - eventually everything should be "maximal" and then loop)
	
	const double diffusionFactor = 0.025;
	const double maxValue = 1;
	controller_ptr->setKused(diffusionFactor);
	controller_ptr->setAused(maxValue);	
	controller_ptr->setDTused(1);	

	auto field_ptr = controller_ptr->getBaseFieldPtr();

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