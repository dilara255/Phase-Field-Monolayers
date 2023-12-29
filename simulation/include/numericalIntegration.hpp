#pragma once

#include <assert.h>

#include "PFM_data.hpp"
#include "derivatives.hpp"
#include "rateOfChangeFunctions.hpp"

//Numerical integration
namespace N_INT { 
	
	enum class rungeKuttaOrder { TWO, FOUR };
	constexpr int rkStepsFromOrder(rungeKuttaOrder order) { 
		if(order == rungeKuttaOrder::TWO) { return 2; }
		else if (order == rungeKuttaOrder::FOUR) { return 4; }
		else { assert("BAD RK ORDER"); return 0; }
	}

	//These correspond to tableus where a_ij = 0 for i != j and c_i = a_i
	double rungeKutaKnBcoef(rungeKuttaOrder order, int n);
	double rungeKutaKnAandCcoef(rungeKuttaOrder order, int n);

//Numerical integration of 2D fields
namespace TD {

//Numerical integration of Cahn-Hilliard potential on 2D fields
namespace CH {
	

void ftcsStep(PFM::PeriodicDoublesLattice2D* phiField, 
			  PFM::PeriodicDoublesLattice2D* auxField,
	          const double dt, const double chK, const double chA, 
	          PFM::checkData_t* checks_ptr);

//2 step Predictor Corrector, a RK2, the trapezoidal version of FTCS
void heunStep(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	          PFM::PeriodicDoublesLattice2D* tempKs_ptr,
	          const double dt, const double chK, const double chA, 
	          PFM::checkData_t* checks_ptr);

//So far only supports order 2 or 4 (with a single predefined tableu for each)
//Will do nothing (or fail an assert) if another oreder is passed
//Only supports tableus where a_ij = 0 for i != j and c_i = a_i
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
