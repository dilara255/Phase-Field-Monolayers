#pragma once

#include "PFM_data.hpp"
#include "derivatives.hpp"

#include "fAux/API/miscStdHeaders.h"

namespace INT {
	//"f", the time derivative for numerical integration, for Cahn-Hiliard
	//TODO: that's not really what this is
	inline double chNumericalF(const PFM::neighborhood9_t* neigh_ptr, const double k, const double A) {
		const double phi = neigh_ptr->getCenter();
		const double laplacianPhi = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

		return k*laplacianPhi - A*phi*(1-phi)*(1-2*phi);
		//+lap -> diffusion (- would accentuate differences)
		//-phi*(1-phi)*(1-2*phi): pulls to the attractors (+ would push)
	}

	/*
	//WARNING: this is not right
	//TODO: FIX IF VERLET IS EVER RE-IMPLEMENTED
	//"A(xn)", the second time derivative for numerical integration, for Cahn-Hiliard
	inline double chNumericalA(const PFM::neighborhood9_t* phiNeigh_ptr, const PFM::neighborhood9_t* dPhiNeigh_ptr, 
																					const double k, const double A) {
	
		const double a = 3;
		const double b = std::sqrt(a);
	
		const double phi = phiNeigh_ptr->getCenter();
		const double dPhi = dPhiNeigh_ptr->getCenter();
		const double laplacianDPhi = PFM::laplacian9pointsAroundNeighCenter(dPhiNeigh_ptr);
		
		return k*laplacianDPhi - A*dPhi*((1 - (a+b)*phi)*(1 - (a-b)*phi));
	}*/
}


namespace INT {
namespace PHU {

	//These are rate of change functions which were wrong but possibly fun. Nothing to see here.

	//"f", the time derivative for numerical integration, for Cahn-Hiliard
	inline double chNumericalF(const PFM::neighborhood9_t* neigh_ptr, const double k, const double A) {
		const double phi = neigh_ptr->getCenter();
		const double laplacian = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

		return k*laplacian - A*phi*(1-phi)*(1-2*phi);
	}

	//"A(xn)", the second time derivative for numerical integration, for Cahn-Hiliard
	//Note: changing the signs of the terms gives off interesting effects
	//-k term + A term seems stable (and fun)
			
	inline double chNumericalA(const PFM::neighborhood9_t* phiNeigh_ptr, const PFM::neighborhood9_t* dPhiNeigh_ptr, 
																					const double k, const double A) {
	
		const double a = 3;
		const double b = std::sqrt(a);
	
		const double phi = phiNeigh_ptr->getCenter();
		const double dPhi = dPhiNeigh_ptr->getCenter();
		const double laplacianDPhi = PFM::laplacian9pointsAroundNeighCenter(dPhiNeigh_ptr);
		
		return k*laplacianDPhi - A*dPhi*((1 - (a+b)*phi)*(1 - (a-b)*phi));

		/*More fun stuff :
		const double laplacianPhi = PFM::laplacian9pointsAroundNeighCenter(phiNeigh_ptr);

		return -k*dPhi*laplacianPhi + A*dPhi*((1 - (a+b)*phi)*(1 - (a-b)*phi));

		return -k + A*dPhi*((1 - (a+b)*phi)*(1 - (a-b)*phi));

		return -k*dPhi*laplacianDPhi + A*dPhi*((1 + (a+b)*phi)*(1 + (a-b)*phi));
		*/
	}

	inline double chNumericalAF(const PFM::neighborhood9_t* neigh_ptr, const double k, const double A) {
		const double phi = neigh_ptr->getCenter();
		const double laplacian = PFM::laplacian9pointsAroundNeighCenter(neigh_ptr);

		const double a = 3;
		const double b = std::sqrt(a);

		double dPhi = k*laplacian - A*phi*(1-phi)*(1-2*phi);

		return k*laplacian - A*dPhi*((1 - (a+b)*phi)*(1 - (a-b)*phi));
	}
}
}