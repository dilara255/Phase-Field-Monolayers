#include "numericalIntegration.hpp"
#include "PFM_dataDefaults.hpp"

typedef struct substep_st {
	double timeAdvanced = 0;
	double currentPhi = 0;
	const double initialPhi = 0;
	double changeRejected = 0;
	const PFM::coordinate_t coordinate;
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

void updateChecks(PFM::checkData_t* checks_ptr, const double change) {

	checks_ptr->densityChange += change;
	checks_ptr->absoluteChangeLastStep += std::abs(change);
	checks_ptr->totalAbsoluteChange += std::abs(change);
	checks_ptr->sumOfsquaresOfChangesLastStep += change * change;
}

namespace N_INT { namespace TD { namespace CH { 
	static std::vector<substep_t> g_substepVector; 
} } }

//Note: this does has quite a bit of coupling with ftcsStepWithSubsteps:
//- It has a little bit of extra logic to also be able to be used when substeps are required. No real overhead;
//- substepAuxField should only ever be used if maySubstep == true;
//- this has to update substepAuxField according to the use in the substepping version:
//-- Currently, the substepping version is using the dPHi of each element. Reason is that this way
//   we don't have to track what substepping elements might be neighbors of other elements also substepping
//   and how much time has passed for them. On the other hand, this makes the initial copy slower, since
//   copying the original values can be done faster, in a single pass (left commented out).
void N_INT::TD::CH::ftcsStep(PFM::PeriodicDoublesLattice2D* phiField, 
						   PFM::PeriodicDoublesLattice2D* auxField,
	                       const double dt, const double chK, const double chA, 
	                       PFM::checkData_t* checks_ptr, 
	                       PFM::PeriodicDoublesLattice2D* substepAuxField, bool maySubstep) {
	
	//if(maySubstep) { substepAuxField->mirrorAllDataFrom(phiField); }

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

	double dPhi, delta;

	//Then we go trhough each point calculating the laplacian of the petential
	for (int j = 0; j < height; j++) {
			centerPoint.y = j;
		for (int i = 0; i < width; i++) {
			centerPoint.x = i;
		
			neigh = auxField->getNeighborhood(centerPoint);
		
			dPhi = -PFM::laplacian9pointsAroundNeighCenter(&neigh);
			//-lap: the more a point "needs" to be pulled than their neighbor, the more they are pulled
			//at the same time, the laplacian has sum 0 (on this periodic field), so total phi is conserverd

			delta = dPhi * dt;

			if (!maySubstep) {
				//We update the field and the checks:
				phiField->incrementDataPoint(centerPoint, delta);
				updateChecks(checks_ptr, delta);
			}
			else {
				//in case substepping is enabled, we need to do some extra book keeping:

				if (std::abs(delta) <= PFM::defaultMaxChangePerStep) {
					//We accept the step and update the checks and dPhis:
					substepAuxField->writeDataPoint(centerPoint, dPhi);
					phiField->incrementDataPoint(centerPoint, delta);
					updateChecks(checks_ptr, delta);
				}
				else {
					//The step was too large, so we'll have to sub-step this
					g_substepVector.push_back({0.0, phiField->getDataPoint(centerPoint), 
						                       phiField->getDataPoint(centerPoint), 
						                       delta, centerPoint});
					//And since no change was made to this element, we'll mark its dPhi as 0:
					substepAuxField->writeDataPoint(centerPoint, 0.0);
				}
			}
		}
	}
}

PFM::neighborhood9_t computeFforAneighborhood9(PFM::neighborhood25_t* neigh25_ptr, 
	                                           const double chK, const double chA) {
		
	PFM::neighborhood9_t fNeigh;
	double numF;

	//We want to calculate 9 values "in the center", each using its neighborhood:
	for (int j = 1; j < 4; j++) { //1, 2, 3: the center
		for (int i = 1; i < 4; i++) {
			
			//But to calculate each value, we actually need the neighborhood around the corresponding point:
			PFM::neighborhood9_t auxNeigh;

			for (int l = -1; l < 2; l++) { //-1, 0, 1: around the center
				for (int k = -1; k < 2; k++) {
					
					const int row = j + l;
					const int collum = i + k;
					
					//also, k and l need to be "recentered" by adding 1:
					auxNeigh.setElement(k+1, l+1, neigh25_ptr->getElement(collum, row));
				}
			}
			
			numF = N_INT::chNumericalF(&auxNeigh, chK, chA);

			const int finalRow = j - 1;
			const int finalCollum = i - 1;

			fNeigh.setElement(finalCollum, finalRow, numF);
		}
	}

	return fNeigh;
}

//TODO: the vector is sorted by coordinate, so we could do a much smarter search. Would that be better?
bool isThisCoordinateOnThisSubstepVector(const PFM::coordinate_t coordinate,
					                     std::vector<substep_t>* const substepVec_ptr) {

	bool found = false;
	int size = substepVec_ptr->size();
	int nextElement = 0;

	while ((!found) && (nextElement < size)) {
		
		found = (substepVec_ptr->at(nextElement).coordinate == coordinate);
		nextElement++;
	}

	return found;
}

//WARNING: the usage of substepAuxField here has to be "synched" with its preparation in ftcsStep.
//TODO: review / test other strategies to improve performance, ie:
//- Could use smaller data structures to hold only the needed information, instead of the complete fields;
//- Receiving the initial phis would simplify the first substep, which is potentially a large chunk of the work;
//- Avoiding re-populating the entire neighborhoods every step;
//- Avoiding multiple computations of the potentials for neighbors of multiple substepping elements;
//- Improvements to make rebalancing easier (and better, maybe).
void N_INT::TD::CH::ftcsStepWithSubsteps(PFM::PeriodicDoublesLattice2D* phiField, 
						   PFM::PeriodicDoublesLattice2D* auxField,
	                       const double dt, const double chK, const double chA, 
	                       PFM::checkData_t* checks_ptr,
	                       PFM::PeriodicDoublesLattice2D* substepAuxField) {
	
	assert(dt > 0);

	//First try the normal step:
	ftcsStep(phiField, auxField, dt, chK, chA, checks_ptr, substepAuxField, true);

	//And then deal with substeps for the elements in need of it.
	
	//note: all times will be treated as if the start of the step was at time 0, since only differences matter
	double smallestTimeAdvanced = 0;
	double lastSmallestTimeAdvance = 2*dt; //large value so the first pass works
	double dtSub = dt/2;
	int extraSubsteps = 0;
	int substepLoops = 0;

	PFM::neighborhood25_t phiNeigh;
	PFM::neighborhood25_t dPhiNeigh;
	PFM::neighborhood9_t fNeigh;

	const size_t totalElements = g_substepVector.size();

	
	while (smallestTimeAdvanced < dt) {
		
		substepLoops++;

		for (size_t i = 0; i < totalElements; i++) {
			
			//At each "pass" here, we'll only actually deal with the least advanced elements, and ignore the others.
			//An alternative to stepping through all elements would be to sort the vector first.
			//TODO: test wether that may actually be faster in realistic use cases.

			//So, in case this element is lagging behind, we try to update it:
			substep_t* const elementToSubstep_ptr = &g_substepVector.at(i);
			const double timeAdvanced = elementToSubstep_ptr->timeAdvanced;
			const double timeToRewind = dt - smallestTimeAdvanced;

			if (timeAdvanced == smallestTimeAdvanced) {
				
				extraSubsteps++;
				const PFM::coordinate_t centerPoint = elementToSubstep_ptr->coordinate;
				
				//auxField has the dPhi for each element and phiField has the last values.
				//Note that in general we do need to update the phiNeigh each step,
				//because of possible interactions with other substepping elements.
				//Also, for elements being substepped, phiField starts with their value at the start of the step.
				//Later, it will hold projections for the final value. This is so we can deal with substepping
				//neighbors the same way we deal with all neighbors. This also works in the first step because
				//both dPhi and time for elements to be substepped start at zero.

				phiNeigh = phiField->getNeighborhood25(centerPoint);
				dPhiNeigh = substepAuxField->getNeighborhood25(centerPoint);

				//We use a linear interpolation to "rewind" the values of the neighborhood to "timeAdvanced":
				dPhiNeigh *= timeToRewind;
				phiNeigh -= dPhiNeigh;

				//Later, we'll store extrapolated intermediate results in the phiField (to help other substeps),
				//so here we need to revert to the actual current value:
				phiNeigh.setCenter(elementToSubstep_ptr->currentPhi);

				//Now we can calculate the new dPhi:
				fNeigh = computeFforAneighborhood9(&phiNeigh, chK, chA);
				const double dPhi = -PFM::laplacian9pointsAroundNeighCenter(&fNeigh);
				double change = dPhi*dtSub;
				
				//And check wether the size of the change was acceptable:
				if (std::abs(change) <= PFM::defaultMaxChangePerStep) {
					
					//Accept the change amd update auxiliar data

					elementToSubstep_ptr->currentPhi += change;
					elementToSubstep_ptr->timeAdvanced += dtSub;
				}

				//If the change was too large, just ignore it (and the next try will be with half the dt)
			}
		}

		//We need to update the smallest value as well as the auxiliar field with the new data:
		smallestTimeAdvanced = 2*dt; //just a large value
		for (size_t i = 0; i < totalElements; i++) {

			auto elementToSubstep_ptr = &g_substepVector.at(i);

			//Update smallest:
			if (g_substepVector.at(i).timeAdvanced < smallestTimeAdvanced) {			
				smallestTimeAdvanced = elementToSubstep_ptr->timeAdvanced;
			}

			//For elements which have stepped, update aux data (even if they haven't stepped this pass)
			//TODO: maybe only update the relevant ones somehow?
			if (elementToSubstep_ptr->timeAdvanced > 0) {
				double effectiveDphi = (elementToSubstep_ptr->currentPhi - elementToSubstep_ptr->initialPhi) 
					                   / elementToSubstep_ptr->timeAdvanced;
				double projectedFinalChange = effectiveDphi * dt;

				phiField->writeDataPoint(elementToSubstep_ptr->coordinate, 
					                     elementToSubstep_ptr->initialPhi + projectedFinalChange);					
				substepAuxField->writeDataPoint(elementToSubstep_ptr->coordinate, effectiveDphi);
			}					
		}

		if (smallestTimeAdvanced == lastSmallestTimeAdvance) {
			//the step was too large for some elements, so:
			dtSub /= 2;
		}
		else { dtSub *= 2; } //the last step size was good, so let's try to pick up speed
		
		lastSmallestTimeAdvance = smallestTimeAdvanced;

		if (dtSub > (dt - smallestTimeAdvanced)) { dtSub = dt - smallestTimeAdvanced; } //prevent over-stepping		
	}

	//The original step had the neighbors of the substepped values change their values based on the
	//initial accessment of the elements which were now substepped. 
	//This new change breaks the simmetry of the laplacian.
	//To restore that we need to add to each neighbor of each element the difference 
	//between the new and the old change. Note, however, that if the neighbor was substepped, than
	//it actually never received the initial change. So we'll need to know which elements are these.
	
	//So this is how we restore that balance (and update the checks):
	for (size_t i = 0; i < totalElements; i++) {
		const double changeAccepted = g_substepVector.at(i).currentPhi - g_substepVector.at(i).initialPhi;
		const double changeRejected = g_substepVector.at(i).changeRejected;
		
		static const double proportionSide = 1/6.0;
		static const double proportionCorner = 1/12.0;
		double proportion;

		double actualDifferenceForNeighbors = 0;
		const int x0 = g_substepVector.at(i).coordinate.x;
		const int y0 = g_substepVector.at(i).coordinate.y;

		for (int j = -1; j < 2; j++) {
			for (int i = -1; i < 2; i++) {
				if (i != 0 || j != 0) { //only neighbors

					const int x = x0+i;
					const int y = y0+j;

					const bool hasSubstepped = isThisCoordinateOnThisSubstepVector({x, y}, &g_substepVector);
					const double differenceOfChanges = changeAccepted - ( (!hasSubstepped) * changeRejected);
					
					//(-1,-1) is a corner, and then it alternates, so:
					if( (i+j)%2 == 0) { proportion = proportionCorner; } 
					else { proportion = proportionSide; };

					const double actualChange = differenceOfChanges * proportion;

					actualDifferenceForNeighbors += actualChange;

					phiField->incrementDataPoint({x, y}, actualChange);
				}
			}
		}

		//For the checks, we need to account both for the final delta phi and for the effect on the neighbors.
		const double deltaPhi = g_substepVector.at(i).currentPhi - g_substepVector.at(i).initialPhi;
		const double effectiveNeighborChange = actualDifferenceForNeighbors;
		updateChecks(checks_ptr, deltaPhi + effectiveNeighborChange);
		//note that this "final" dPhi may actually be final, since neighbors may change it here as well.
	}
	
	//TODO: when restarting the simulation via GUI, these vectors will be as large as the largest ever needed
	//This might mean unecessary memory usage: profile and decide wether that's okay.
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