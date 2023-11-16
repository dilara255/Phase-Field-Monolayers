#pragma once

#include "PFM_data.hpp"
#include "derivatives.hpp"

#include "fAux/API/miscStdHeaders.h"

namespace INT {
	//"f", the time derivative for numerical integration, for Cahn-Hiliard
	inline double chNumericalF(const PFM::neighborhood9_t* neigh_ptr, const double k, const double A) {
		const double phi = neigh_ptr->getCenter();
		const double laplacian = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

		return k*laplacian - A*phi*(1-phi)*(1-2*phi);
	}

	//"A(xn)", the second time derivative for numerical integration, for Cahn-Hiliard
	inline double chNumericalA(const PFM::neighborhood9_t* phiNeigh_ptr, const PFM::neighborhood9_t* dPhiNeigh_ptr, 
																					const double k, const double A) {
	
		const static double a = 3 + std::sqrt(3);
		const static double b = 3 - std::sqrt(3);
	
		const double phi = phiNeigh_ptr->getCenter();
		const double dPhi = dPhiNeigh_ptr->getCenter();
		const double laplacianDPhi = PFM::laplacian9pointsAroundNeighCenter(dPhiNeigh_ptr);

		return k*laplacianDPhi - A*dPhi*( (1 - a*phi)*(1 - b*phi) + 1);
	}
}