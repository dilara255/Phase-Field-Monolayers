#pragma once

#include <assert.h>

#include "PFM_data.hpp"
#include "derivatives.hpp"
#include "rateOfChangeFunctions.hpp"

namespace INT { 
	
	enum class rungeKuttaOrder { TWO, FOUR };
	inline int rkStepsFromOrder(rungeKuttaOrder order) { 
		if(order == rungeKuttaOrder::TWO) { return 2; }
		else if (order == rungeKuttaOrder::FOUR) { return 4; }
		else { assert("BAD RK ORDER"); return 0; }
	}

	double rungeKutaKnFinalCoef(rungeKuttaOrder order, int n);

	double rungeKutaKnIntermediateCoef(rungeKuttaOrder order, int n);

namespace TD {

void explicitEulerCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                          const double dt, const double chK, const double chA, 
	                          PFM::checkData_t* checks_ptr);

//For now, runs for a set amount of substeps
//TODO: actually test convergence
void implicitEulerCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                          PFM::PeriodicDoublesLattice2D* baseField_ptr, const double dt,
	                          const double chK, const double chA, PFM::checkData_t* checks_ptr);

void heunCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	                 PFM::PeriodicDoublesLattice2D* tempKs_ptr,
	                 const double dt, const double chK, const double chA, 
	                 PFM::checkData_t* checks_ptr);

void verletCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA, 
	                               PFM::checkData_t* checks_ptr);

//So far only supports 2 or 4 steps. Will do nothing (or fail an assert) if another amount is 
//TODO: add more coeficient lists and order enumerations
//TODO: eventually, maybe, generalize or something, idk
void rungeKuttaCahnHiliard(rungeKuttaOrder order, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
						  PFM::PeriodicDoublesLattice2D* field_ptr,
	                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                      const double dt, const double chK, const double chA, PFM::checkData_t* checks_ptr);

}
}
