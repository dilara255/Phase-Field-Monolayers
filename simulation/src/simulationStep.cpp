#include "simulControl.hpp"
#include "derivatives.hpp"
#include <assert.h>

//TODO: make an initial softening (that works)
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

void PFM::singleLayerSim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//Test simulation: the pixels initialized as non-zero should "diffuse" up and to the right
	//(not really diffuse, more like reinforce - eventually everything should be "maximal" and then loop)
	
	const double k = 0.5;
	const double A = 0.01;
	const double dt = 0.05;
	const double initialCellDiameterDensity = 1; //always a bit less because integers. Also, seeds not packed
	//TODO: on the seeding function, define spacing (ie, option to use a spation configuration for best packing)

	auto field_ptr = controller_ptr->getFieldPtr();

	const int width = (int)field_ptr->getFieldDimensions().width;
	const int height = (int)field_ptr->getFieldDimensions().height;

	const double initialCellRadius = initialCellDiameterDensity * std::min(width, height)
		                             / (2 * std::sqrt(2) * std::sqrt((float)controller_ptr->getNumberCells()));

	expandCells(field_ptr, initialCellRadius, controller_ptr->getLastCellSeedValue(), true);
	
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

	auto field_ptr = controller_ptr->getFieldPtr();

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