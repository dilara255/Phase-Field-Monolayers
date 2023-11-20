#include "numericalIntegration.hpp"

//TODO: Implement : )

double INT::rungeKutaKnIntermediateCoef(rungeKuttaOrder order, int n) {
	if(n <= 0) { assert(false); return 0.0; }

	double order2[2] = { 0.5, 0 };
	double order4[4] = { 0.5, 0.5, 1, 0 };

	if (order == rungeKuttaOrder::TWO) {
		if(n > 2) { assert(false); return 0.0; }

		return order2[n-1];
	}
	else if (order == rungeKuttaOrder::FOUR) {
		if(n > 4) { assert(false); return 0.0; }

		return order4[n-1];
	}
	else { assert(false); return 0.0; }
}

double INT::rungeKutaKnFinalCoef(rungeKuttaOrder order, int n) {
	if(n <= 0) { assert(false); return 0.0; }

	double order2[2] = { 0, 1 };
	double order4[4] = { 1.0/6, 1.0/3, 1.0/3, 1.0/6 };

	if (order == rungeKuttaOrder::TWO) {
		if(n > 2) { assert(false); return 0.0; }

		return order2[n-1];
	}
	else if (order == rungeKuttaOrder::FOUR) {
		if(n > 4) { assert(false); return 0.0; }

		return order4[n-1];
	}
	else { assert(false); return 0.0; }
}

//Really this is FTCS
void INT::TD::explicitEulerCahnHiliard(PFM::PeriodicDoublesLattice2D* phiField, 
									   PFM::PeriodicDoublesLattice2D* auxField,
	                                   const double dt, const double chK, const double chA, 
	                                   PFM::checkData_t* checks_ptr) {
	
	int height = phiField->getFieldDimensions().height;
	int width = phiField->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = phiField->getNeighborhood(centerPoint);
		
			auxField->writeDataPoint(centerPoint, chNumericalF(&neigh, chK, chA));
		}
	}

	double dPhi;

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = auxField->getNeighborhood(centerPoint);
		
			dPhi = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
			//-lap: the more a point "needs" to be pulled than their neighbor, the more they are pulled
			//at the same time, the laplacian has sum 0 (on this periodic field), so total phi is conserverd

			phiField->incrementDataPoint(centerPoint, dt * dPhi);
			
			checks_ptr->densityChange += dPhi * dt;
			checks_ptr->absoluteChange += std::abs(dPhi);
		}
	}
}

//TODO: FIX ALL OF THE BELLOW:
/*
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

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	int height = currentStepField_ptr->getFieldDimensions().height;
	int width = currentStepField_ptr->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t phiNeigh;
	PFM::neighborhood9_t dPhiNeigh;

	double an, phi0, dPhi0, phi, phiPrev; //phiPrev: phi_n-1

	//dPhi0 = f(phi0)
	//an = A(phi0, dPhi0)
	//phi = 2*phi0 - phiPrev + an * dt * dt

	//First we want to get dPhi on all elements, plus store their current phi
	//The current phi will later become the "next step's previous step's phi", after the fields rotate (in the end)
	//Note along the way that we only write to currentStepField_ptr and only read lastStepField_ptr
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			phiNeigh = field_ptr->getNeighborhood(centerPoint);
			phi0 = phiNeigh.getCenter();

			dPhi0 = chNumericalF(&phiNeigh, chK, chA);
			tempKs_ptr->writeDataPoint(centerPoint, dPhi0);

			//Note: currentStepField is written before it's read: value comming in doesn't matter
			//After rotation, this will become lastStepField and hold the phiPrev of the next step
			currentStepField_ptr->writeDataPoint(centerPoint, phi0); 
		}
	}

	//Then we can calculate an and phi in a single loop:
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			phiNeigh = field_ptr->getNeighborhood(centerPoint);
			phi0 = phiNeigh.getCenter();

			dPhiNeigh = tempKs_ptr->getNeighborhood(centerPoint);

			an = INT::chNumericalA(&phiNeigh, &dPhiNeigh, chK, chA);

			phiPrev = lastStepField_ptr->getDataPoint(centerPoint); //from the step before the last one completed

			phi = 2*phi0 - phiPrev + an * dt * dt;

			field_ptr->writeDataPoint(centerPoint, phi);

			checks_ptr->density += phi;
			checks_ptr->absoluteChange += std::abs( (phi-phi0)/dt );
		}
	}

	//As noted, we only write to currentStepField_ptr and only read lastStepField_ptr
	//The rotation will make lastStepField take the value of currentStepField
	//which is the value on the step before this one
	//this will be the value of phi_n-1, or phiPrev, in the next step
	rotatingField_ptr->rotatePointers();
}

void rungeKutaCahnHiliardFirstStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA) {

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double phi0, kn;

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();

	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
			
			neigh = field_ptr->getNeighborhood(centerPoint);
				
			kn = dt * INT::chNumericalF(&neigh, chK, chA);
			tempKs_ptr->writeDataPoint(centerPoint, coefKnFinal * kn); //also clears old data
				
			phi0 = neigh.getCenter();
			currentStepField_ptr->writeDataPoint(centerPoint, phi0 + coefKnInterm * kn);
		}
	}
	rotatingField_ptr->rotatePointers();

}

void rungeKutaCahnHiliardIntermediateStep(double coefKnFinal, double coefKnInterm,
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
			
				neigh = lastStepField_ptr->getNeighborhood(centerPoint);

				kn = dt * INT::chNumericalF(&neigh, chK, chA);

				tempKs_ptr->incrementDataPoint(centerPoint, coefKnFinal * kn);
				
				phi0 = field_ptr->getDataPoint(centerPoint);
				currentStepField_ptr->writeDataPoint(centerPoint, phi0 + coefKnInterm * kn);
			}
	}
	rotatingField_ptr->rotatePointers();
}

void rungeKutaCahnHiliardFinalStep(double coefKnFinal, int height, int width, 
								   PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
								   PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
								   const double dt, const double chK, const double chA,
								   PFM::checkData_t* checks_ptr) {
	
	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double dPhi, phi0, phi, kn;

	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
			
			neigh = lastStepField_ptr->getNeighborhood(centerPoint);

			kn = dt * INT::chNumericalF(&neigh, chK, chA);
			dPhi = tempKs_ptr->getDataPoint(centerPoint) + coefKnFinal * kn;

			phi0 = field_ptr->getDataPoint(centerPoint);
			phi = phi0 + dPhi;

			field_ptr->writeDataPoint(centerPoint, phi);		

			checks_ptr->density += phi;
			checks_ptr->absoluteChange += std::abs(dPhi/dt);
		}
	}
}

void INT::TD::rungeKuttaCahnHiliard(INT::rungeKuttaOrder order, 
	                                PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								    PFM::PeriodicDoublesLattice2D* field_ptr,
								    PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
								    const double dt, const double chK, const double chA, 
								    PFM::checkData_t* checks_ptr) {

	int steps = INT::rkStepsFromOrder(order);
	if(steps == 0) { assert(false); return; }

	int height = field_ptr->getFieldDimensions().height;
	int width = field_ptr->getFieldDimensions().width;

	double coefKnInterm = rungeKutaKnIntermediateCoef(order, 1);
	double coefKnFinal = rungeKutaKnFinalCoef(order, 1);
	rungeKutaCahnHiliardFirstStep(coefKnFinal, coefKnInterm, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA);
	if (steps > 2) {
		for(int i = 0; i < (steps - 2); i++) {

			coefKnInterm = rungeKutaKnIntermediateCoef(order, i+2);
			coefKnFinal = rungeKutaKnFinalCoef(order, i+2);
			rungeKutaCahnHiliardIntermediateStep(coefKnFinal, coefKnInterm, height, width, 
				                                 rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA);
		}
	}

	coefKnFinal = rungeKutaKnFinalCoef(order, steps);
	rungeKutaCahnHiliardFinalStep(coefKnFinal, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA, checks_ptr);
}
*/