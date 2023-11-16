#include "numericalIntegration.hpp"

//TODO: Implement : )


double INT::rungeKutaKnFinalCoef(int steps, int n) {
	return 0;
}

double INT::rungeKutaKnIntermediateCoef(int steps, int n) {
	return 0;
}

void INT::TD::explicitEulerCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                                   const double dt, const double chK, const double chA, 
	                                   PFM::checkData_t* checks_ptr) {
	
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
	rotatingField_ptr->rotatePointers();
}

//TODO: pass maxSteps plus maxDif, run while(dif > maxDif && step < maxSteps)
void INT::TD::implicitEulerCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                                   PFM::PeriodicDoublesLattice2D* baseField_ptr, const double dt,
	                                   const double chK, const double chA, PFM::checkData_t* checks_ptr) {

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

void INT::TD::heunCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
							  PFM::PeriodicDoublesLattice2D* tempKs_ptr,
							  const double dt, const double chK, const double chA,
							  PFM::checkData_t* checks_ptr) {

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	int height = currentStepField_ptr->getFieldDimensions().height;
	int width = currentStepField_ptr->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double dPhi0;
	double phi0;
	double dPhi;
	double phi;

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = lastStepField_ptr->getNeighborhood(centerPoint);
			phi0 = neigh.getCenter();
		
			dPhi0 = chNumericalF(&neigh, chK, chA);
			tempKs_ptr->writeDataPoint(centerPoint, dPhi0);

			currentStepField_ptr->writeDataPoint(centerPoint, phi0 + dt * dPhi0);
		}
	}

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			phi0 = lastStepField_ptr->getDataPoint(centerPoint);
			neigh = currentStepField_ptr->getNeighborhood(centerPoint);
		
			dPhi = chNumericalF(&neigh, chK, chA);
			dPhi0 = tempKs_ptr->getDataPoint(centerPoint);

			phi = phi0 + 0.5 * dt * (dPhi + dPhi0);
			currentStepField_ptr->writeDataPoint(centerPoint, phi);

			checks_ptr->density += phi;
			checks_ptr->absoluteChange += std::abs(0.5*(dPhi + dPhi0));
		}
	}
	
	rotatingField_ptr->rotatePointers();
}

void INT::TD::verletCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA, 
	                               PFM::checkData_t* checks_ptr) {

}

void INT::TD::rungeKutaCahnHiliardFirstStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA) {

}

void INT::TD::rungeKutaCahnHiliardIntermediateStep(int RKstep, double coefKnFinal, double coefKnInterm,
	                                      int height, int width, 
	                                      PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
										  PFM::PeriodicDoublesLattice2D* field_ptr,
	                                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                                      const double dt, const double chK, const double chA) {

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

void INT::TD::rungeKutaCahnHiliardFinalStep(double coefKnFinal, double coefKnInterm, int height, int width, 
										    PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
										    PFM::PeriodicDoublesLattice2D* field_ptr,
										    PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
										    const double dt, const double chK, const double chA,
										    PFM::checkData_t* checks_ptr) {

}

void INT::TD::rungeKutaCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
								   PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
								   const double dt, const double chK, const double chA, 
								   PFM::checkData_t* checks_ptr) {

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