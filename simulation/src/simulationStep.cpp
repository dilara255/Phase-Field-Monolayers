#include "simulControl.hpp"

void PFM::singleLayerSim_fn(SimulationControl* controller_ptr, int* stepCount_ptr, bool* isRunning_ptr) {
	//Test simulation: the pixels initialized as non-zero should "diffuse" up and to the right
	//(not really diffuse, more like reinforce - eventually everything should be "maximal" and then loop)
	
	const double diffusionFactor = 0.025;
	const double maxValue = 1;

	auto field_ptr = controller_ptr->getFieldPtr();

	const int width = (int)field_ptr->getFieldDimensions().width;
	const int height = (int)field_ptr->getFieldDimensions().height;
	
	double value, valueBellow, valueToTheLeft;

	while(!controller_ptr->checkIfShouldStop()) {
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {

				value = field_ptr->getDataPoint({i,j});
				valueBellow = field_ptr->getDataPoint({i,j+1});
				valueToTheLeft = field_ptr->getDataPoint({i-1,j});

				field_ptr->writeDataPoint({i,j+1}, valueBellow + (diffusionFactor * value) );
				field_ptr->writeDataPoint({i-1,j}, valueToTheLeft + (diffusionFactor * value) );

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