#include "numericalIntegration.hpp"
#include "PFM_defaults.hpp"

typedef struct substep_st {
	bool done;
	PFM::coordinate_t position;
	double dtNextSubstep;
	double timeRemaining;
} substep_t;

double N_INT::rungeKutaKnAandCcoef(rungeKuttaOrder order, int n) {
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

double N_INT::rungeKutaKnBcoef(rungeKuttaOrder order, int n) {
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

void updateChecks(PFM::checkData_t* checks_ptr, const double dt, const double dPhi) {

	checks_ptr->densityChange += dPhi * dt;
	checks_ptr->absoluteChange += std::abs(dPhi) * dt;
	checks_ptr->sumOfsquaresOfChanges += (dPhi * dt) * (dPhi * dt);
}

namespace N_INT { namespace TD { namespace CH { static std::vector<substep_t> g_substepVector; } } }

//Has a little bit of extra logic to also be able to be used when substeps are required. No real overhead.
void N_INT::TD::CH::ftcsStep(PFM::PeriodicDoublesLattice2D* phiField, 
						   PFM::PeriodicDoublesLattice2D* auxField,
	                       const double dt, const double chK, const double chA, 
	                       PFM::checkData_t* checks_ptr, bool maySubstep) {
	
	int height = phiField->getFieldDimensions().height;
	int width = phiField->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	
	//auxField holds the data of the CH potential for each site
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = phiField->getNeighborhood(centerPoint);
		
			auxField->writeDataPoint(centerPoint, chNumericalF(&neigh, chK, chA));
		}
	}

	double dPhi;

	//Then we go trhough each point calculating the laplacian of the petential
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = auxField->getNeighborhood(centerPoint);
		
			dPhi = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
			//-lap: the more a point "needs" to be pulled than their neighbor, the more they are pulled
			//at the same time, the laplacian has sum 0 (on this periodic field), so total phi is conserverd

			if (!maySubstep || std::abs(dPhi) <= PFM::maxChangePerStep) {
				//We accept the step and update the checks
				phiField->incrementDataPoint(centerPoint, dt * dPhi);
				updateChecks(checks_ptr, dt, dPhi);
			}		
			else {
				//The step was too large, so we'll have to sub-step this
				g_substepVector.push_back({false, centerPoint, dt, dt});
			}
		}
	}
}

void N_INT::TD::CH::ftcsStepWithSubsteps(PFM::PeriodicDoublesLattice2D* phiField, 
						   PFM::PeriodicDoublesLattice2D* auxField,
	                       const double dt, const double chK, const double chA, 
	                       PFM::checkData_t* checks_ptr) {
	
	//First try the normal step:
	ftcsStep(phiField, auxField, dt, chK, chA, checks_ptr, true);

	//And then deal with substeps
	int elementsPending = g_substepVector.size();
	double timestep;
	double dPhi;
	substep_t elementToSubstep;
	double extraSubsteps = 0;
	PFM::coordinate_t centerPoint;
	
	while (elementsPending > 0) {
		
		for (int i = 0; i < elementsPending; i++) {
		
			elementToSubstep = g_substepVector.at(i);

			if (!elementToSubstep.done) {
				//First we check that we won't overstep:
				elementToSubstep.dtNextSubstep /= 2.0;
				timestep = std::min(elementToSubstep.dtNextSubstep, elementToSubstep.timeRemaining);
				
				//Then we try a new substep:
				extraSubsteps++;

				//auxField has the original chNumericalF for all elements
				//phiField has the final value for all elements (so far, or something)
				dPhi = 1.0; //hopefully not : p

				//And check wether the size of the change was acceptable:
				if (std::abs(dPhi) <= PFM::maxChangePerStep) {
					//Accept changes and update checks
					phiField->incrementDataPoint(centerPoint, dt * dPhi);
					updateChecks(checks_ptr, dt, dPhi);

					//Prepare for the next substep:
					elementToSubstep.timeRemaining -= timestep;
					elementToSubstep.dtNextSubstep *= 2.0;
					//But also check for the end of the step:
					if (elementToSubstep.timeRemaining <= 0) {
						elementToSubstep.done = true;
						elementsPending--;
					}					
				}

				//If the change was too large, just ignore it and the next try will be with half the dt
			}
		}
	}

	//TODO: when restarting the simulation via GUI, this vector will be as large as the largest ever needed
	//This might mean unecessary memory usage. Profile and decide wether that's okay.
	g_substepVector.clear();
}

//Essentially a RK2, or a 2 step predictor-corrector, or the trapezoidal rule applied to FTCS
void N_INT::TD::CH::heunStep(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
						   PFM::PeriodicDoublesLattice2D* tempKs_ptr,
						   const double dt, const double chK, const double chA,
						   PFM::checkData_t* checks_ptr) {

	// {1, 0.5, 0.5} or {2/3, 1/4, 3/4} are typical. b1 + b2 should = 1
	const double a1 = 1; 
	const double b1 = 0.5; 
	const double b2 = 0.5;

	auto currentStepField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto lastStepField_ptr = rotatingField_ptr->getPointerToLast();

	int height = currentStepField_ptr->getFieldDimensions().height;
	int width = currentStepField_ptr->getFieldDimensions().width;

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double dPhiIntermediates; //dPhi will be calculated for two points
	double phi0;
	double phi;
	double phiIntermediate;

	//last: last value of phi, phi_n
	//current: effectively junk
	//tempKs_ptr: effectively junk
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = lastStepField_ptr->getNeighborhood(centerPoint);
		
			tempKs_ptr->writeDataPoint(centerPoint, chNumericalF(&neigh, chK, chA));
		}
	}

	//last: last value of phi, phi_n
	//current: effectively junk
	//tempKs_ptr: the inner part of the intermediate dPhi (outer laplacian still to be apllied)
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = tempKs_ptr->getNeighborhood(centerPoint);
			dPhiIntermediates = -PFM::laplacian9pointsAroundNeighCenter(&neigh);

			phi0 = lastStepField_ptr->getDataPoint(centerPoint);
			phiIntermediate = phi0 + a1 * dt * dPhiIntermediates;

			//currentStepField_ptr wiil hold the intermediate value of phi for now
			currentStepField_ptr->writeDataPoint(centerPoint, phiIntermediate);

			//We can "half-step" the original phi value now (to avoid am extra loop):
			lastStepField_ptr->incrementDataPoint(centerPoint, b1 * dt * dPhiIntermediates);

			//And also the checks:
			checks_ptr->densityChange += b1 * dt * dPhiIntermediates;
			checks_ptr->absoluteChange += std::abs(b1 * dPhiIntermediates) * dt; //TODO: is this right? are Kn always monotone?
			checks_ptr->sumOfsquaresOfChanges += (b1 * dPhiIntermediates * dt) * (b1 * dPhiIntermediates * dt); //See here as well
		}
	}

	//last: phi_n + 0.5 * dPhiIntermediate
	//current: phiIntermediate
	//tempKs_ptr: effectively back to junk
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = currentStepField_ptr->getNeighborhood(centerPoint);
		
			tempKs_ptr->writeDataPoint(centerPoint, chNumericalF(&neigh, chK, chA));
		}
	}

	//last: phi_n + 0.5 * dPhiIntermediate
	//current: effectively back to junk
	//tempKs_ptr: the inner part of the SECOND intermediate dPhi (outer laplacian still to be apllied)
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = tempKs_ptr->getNeighborhood(centerPoint);
			dPhiIntermediates = -PFM::laplacian9pointsAroundNeighCenter(&neigh);

			phi = lastStepField_ptr->getDataPoint(centerPoint) + b2 * dt * dPhiIntermediates;
			
			//CurrentStepField_ptr wiil now hold the final value of phi, as expected
			currentStepField_ptr->writeDataPoint(centerPoint, phi);

			//The checks can now be finalized:
			checks_ptr->densityChange += b2 * dt * dPhiIntermediates;
			checks_ptr->absoluteChange += std::abs(b2 * dPhiIntermediates) * dt;
			checks_ptr->sumOfsquaresOfChanges += (b2 * dPhiIntermediates * dt) * (b2 * dPhiIntermediates * dt);
		}
	}
	
	rotatingField_ptr->rotatePointers();
}

//The RK methods implemented bellow pressupose tableaus where a_ij = 0 for i != j and c_i = a_i
//TODO: check whether that's enough

void rungeKutaCahnHiliardFirstStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA,
								   PFM::checkData_t* checks_ptr) {

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double phi0, k1, dPhiIntermediate;

	auto newPhiField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto tempForDphiField_ptr = rotatingField_ptr->getPointerToLast();
	//field_ptr: last value of phi, phi_n

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = field_ptr->getNeighborhood(centerPoint);
		
			tempForDphiField_ptr->writeDataPoint(centerPoint, N_INT::chNumericalF(&neigh, chK, chA));
		}
	}

	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
			
			neigh = tempForDphiField_ptr->getNeighborhood(centerPoint);
				
			k1 = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
			dPhiIntermediate = dt * k1;

			phi0 = field_ptr->getDataPoint(centerPoint);

			//Will hold phi_n + a1 * h * k1
			tempKs_ptr->writeDataPoint(centerPoint, phi0 + coefKnInterm * dPhiIntermediate); 

			//currentStepField_ptr wiil hold the intermediate steps for the final value of phi
			//Since the calculation of phi is linear on the intermediate Ks, we can do it in parts
			//Since we're using tableus with a_ij = 0 for i =! j, that means we just need to keep one Kn per step
			newPhiField_ptr->writeDataPoint(centerPoint, phi0 + coefKnFinal * dPhiIntermediate);

			//And also the checks:
			checks_ptr->densityChange += coefKnFinal * dPhiIntermediate;
			checks_ptr->absoluteChange += std::abs(coefKnFinal * k1) * dt; //TODO: is this right? are Kn always monotone?
			checks_ptr->sumOfsquaresOfChanges += (coefKnFinal * k1 * dt) * (coefKnFinal * k1 * dt);
		}
	}
}

void rungeKutaCahnHiliardIntermediateStep(double coefKnFinal, double coefKnInterm,
	                                      int height, int width, 
	                                      PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
										  PFM::PeriodicDoublesLattice2D* field_ptr,
	                                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                                      const double dt, const double chK, const double chA,
								          PFM::checkData_t* checks_ptr) {

	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double phi0, kn, dPhiIntermediate;

	auto newPhiField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto tempForDphiField_ptr = rotatingField_ptr->getPointerToLast();
	//field_ptr: last value of phi, phi_n

	//tempKs: holds the effective phi to be used when calculating the next K
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = tempKs_ptr->getNeighborhood(centerPoint);
		
			tempForDphiField_ptr->writeDataPoint(centerPoint, N_INT::chNumericalF(&neigh, chK, chA));
		}
	}

	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

			for (int i = 0; i < width; i++) {
				centerPoint.x = i;
			
				neigh = tempForDphiField_ptr->getNeighborhood(centerPoint);
				
				kn = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
				dPhiIntermediate = dt * kn;

				phi0 = field_ptr->getDataPoint(centerPoint);

				//Will hold phi_n + a1 * h * kn
				tempKs_ptr->writeDataPoint(centerPoint, phi0 + coefKnInterm * dPhiIntermediate); 

				//currentStepField_ptr wiil hold the intermediate steps for the final value of phi
				//Since the calculation of phi is linear on the intermediate Ks, we can do it in parts
				//Since we're using tableus with a_ij = 0 for i =! j, that means we just need to keep one Kn per step
				newPhiField_ptr->incrementDataPoint(centerPoint, coefKnFinal * dPhiIntermediate);

				//And also the checks:
				checks_ptr->densityChange += coefKnFinal * dPhiIntermediate;
				checks_ptr->absoluteChange += std::abs(coefKnFinal * kn) * dt; //TODO: is this right? are Kn always monotone?
				checks_ptr->sumOfsquaresOfChanges += (coefKnFinal * kn * dt) * (coefKnFinal * kn * dt);
			}
	}
}

void rungeKutaCahnHiliardFinalStep(double coefKnFinal, int height, int width, 
								   PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
								   PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
								   const double dt, const double chK, const double chA,
								   PFM::checkData_t* checks_ptr) {
	
	PFM::coordinate_t centerPoint;
	PFM::neighborhood9_t neigh;
	double kn, dPhiIntermediate;

	auto newPhiField_ptr = rotatingField_ptr->getPointerToCurrent();
	auto tempForDphiField_ptr = rotatingField_ptr->getPointerToLast();
	//field_ptr: last value of phi, phi_n

	//tempKs: holds the effective phi to be used when calculating the next K
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = tempKs_ptr->getNeighborhood(centerPoint);
		
			tempForDphiField_ptr->writeDataPoint(centerPoint, N_INT::chNumericalF(&neigh, chK, chA)); //also clears old data
		}
	}

	for (int j = 0; j < height; j++) {
		centerPoint.y = j;

		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
			
			neigh = tempForDphiField_ptr->getNeighborhood(centerPoint);

			kn = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
			dPhiIntermediate = dt * kn;
		
			//currentStepField_ptr wiil now hold the final new value of phi
			newPhiField_ptr->incrementDataPoint(centerPoint, coefKnFinal * dPhiIntermediate);

			field_ptr->writeDataPoint(centerPoint, newPhiField_ptr->getDataPoint(centerPoint));

			//And also the checks:
			checks_ptr->densityChange += coefKnFinal * dPhiIntermediate;
			checks_ptr->absoluteChange += std::abs(coefKnFinal * kn) * dt; //TODO: is this right? are Kn always monotone?
			checks_ptr->sumOfsquaresOfChanges += (coefKnFinal * kn * dt) * (coefKnFinal * kn * dt);
		}
	}
}

void N_INT::TD::CH::rungeKuttaStep(N_INT::rungeKuttaOrder order, 
	                             PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								 PFM::PeriodicDoublesLattice2D* field_ptr,
								 PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
								 const double dt, const double chK, const double chA, 
								 PFM::checkData_t* checks_ptr) {

	int steps = N_INT::rkStepsFromOrder(order);
	if(steps == 0) { assert(false); return; }

	int height = field_ptr->getFieldDimensions().height;
	int width = field_ptr->getFieldDimensions().width;

	double coefKnInterm = rungeKutaKnAandCcoef(order, 1);
	double coefKnFinal = rungeKutaKnBcoef(order, 1);
	rungeKutaCahnHiliardFirstStep(coefKnFinal, coefKnInterm, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA, checks_ptr);
	if (steps > 2) {
		for(int i = 0; i < (steps - 2); i++) {

			coefKnInterm = rungeKutaKnAandCcoef(order, i+2);
			coefKnFinal = rungeKutaKnBcoef(order, i+2);
			rungeKutaCahnHiliardIntermediateStep(coefKnFinal, coefKnInterm, height, width, 
				                                 rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA, checks_ptr);
		}
	}

	coefKnFinal = rungeKutaKnBcoef(order, steps);
	rungeKutaCahnHiliardFinalStep(coefKnFinal, height, width, 
		                          rotatingField_ptr, field_ptr, tempKs_ptr, dt, chK, chA, checks_ptr);
}

/*TODO: Decide wether or not to actually implement this one
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
			checks_ptr->densityChange += phi-phi0;
			checks_ptr->absoluteChange += std::abs( phi-phi0 );
			checks_ptr->sumOfsquaresOfChanges += (phi-phi0) * (phi-phi0);
		}
	}

	//As noted, we only write to currentStepField_ptr and only read lastStepField_ptr
	//The rotation will make lastStepField take the value of currentStepField
	//which is the value on the step before this one
	//this will be the value of phi_n-1, or phiPrev, in the next step
	rotatingField_ptr->rotatePointers();
}
*/