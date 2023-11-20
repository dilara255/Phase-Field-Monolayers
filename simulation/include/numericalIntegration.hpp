#pragma once

#include <assert.h>

#include "PFM_data.hpp"
#include "derivatives.hpp"
#include "rateOfChangeFunctions.hpp"

//Numerical integration
namespace INT { 
	
	enum class rungeKuttaOrder { TWO, FOUR };
	inline int rkStepsFromOrder(rungeKuttaOrder order) { 
		if(order == rungeKuttaOrder::TWO) { return 2; }
		else if (order == rungeKuttaOrder::FOUR) { return 4; }
		else { assert("BAD RK ORDER"); return 0; }
	}

	double rungeKutaKnFinalCoef(rungeKuttaOrder order, int n);

	double rungeKutaKnIntermediateCoef(rungeKuttaOrder order, int n);

//Numerical integration of 2D fields
namespace TD {

//Numerical integration of Cahn-Hilliard potential on 2D fields
namespace CH {
	

void fctsStep(PFM::PeriodicDoublesLattice2D* phiField, 
			  PFM::PeriodicDoublesLattice2D* auxField,
	          const double dt, const double chK, const double chA, 
	          PFM::checkData_t* checks_ptr);

void heunStep(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	          PFM::PeriodicDoublesLattice2D* tempKs_ptr,
	          const double dt, const double chK, const double chA, 
	          PFM::checkData_t* checks_ptr);

//So far only supports order 2 or 4. Will do nothing (or fail an assert) if another oreder is passed
//TODO: add more coeficient lists and order enumerations
//TODO: eventually, maybe, generalize or something, idk
void rungeKuttaStep(rungeKuttaOrder order, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
				    PFM::PeriodicDoublesLattice2D* field_ptr,
	                PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                const double dt, const double chK, const double chA, PFM::checkData_t* checks_ptr);

/* TODO: either reimplement or clear this
void verletCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               const double dt, const double chK, const double chA, 
	                               PFM::checkData_t* checks_ptr);
*/
}
}
}
