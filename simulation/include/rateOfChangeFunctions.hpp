#pragma once

#include "PFM_data.hpp"
#include "derivatives.hpp"

//"f", the time derivative for numerical integration, for Cahn-Hiliard
inline double chNumericalF(const PFM::neighborhood9_t* neigh_ptr, double k, double A) {
	double phi = neigh_ptr->getCenter();
	double laplacian = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

	return k*laplacian - A*phi*(1-phi)*(1-2*phi);
}