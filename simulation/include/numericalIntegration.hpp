#pragma once

#include <assert.h>

#include "PFM_data.hpp"
#include "derivatives.hpp"
#include "rateOfChangeFunctions.hpp"

namespace INT { 
	
	double rungeKutaKnFinalCoef(int steps, int n);

	double rungeKutaKnIntermediateCoef(int steps, int n);

namespace TD {

void explicitEulerCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                          double dt, double chK, double chA, PFM::checkData_t* checks_ptr);

//TODO: maxSteps plus maxDif, e fazer rodar até (dif <= maxDif || step >= maxSteps)
void implicitEulerCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr, 
	                                            PFM::PeriodicDoublesLattice2D* baseField_ptr, double dt,
	                                               double chK, double chA, PFM::checkData_t* checks_ptr);

void heunCahnHiliard(PFM::PeriodicDoublesLattice2D* field_ptr,
	                 PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	                 double chK, double chA, double dt, 
	                 PFM::checkData_t* checks_ptr);

void verletCahnHiliard(PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               double chK, double chA, double dt, 
	                               PFM::checkData_t* checks_ptr);

void rungeKutaCahnHiliardFirstStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
								   PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               double chK, double chA, double dt);

void rungeKutaCahnHiliardIntermediateStep(int RKstep, double coefKnFinal, double coefKnInterm,
	                                      int height, int width, 
	                                      PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
										  PFM::PeriodicDoublesLattice2D* field_ptr,
	                                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                                      double chK, double chA, double dt);

void rungeKutaCahnHiliardFinalStep(double coefKnFinal, double coefKnInterm, int height, int width, 
	                               PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
	                               PFM::PeriodicDoublesLattice2D* field_ptr,
	                               PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                               double chK, double chA, double dt, PFM::checkData_t* checks_ptr);

void rungeKutaCahnHiliard(int steps, PFM::CurrentAndLastPerioricDoublesLattice2D* rotatingField_ptr,
						  PFM::PeriodicDoublesLattice2D* field_ptr,
	                      PFM::PeriodicDoublesLattice2D* tempKs_ptr, 
	                      double chK, double chA, double dt, PFM::checkData_t* checks_ptr);

}
}
