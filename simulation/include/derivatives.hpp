#pragma once

#include "PFM_data.hpp"

namespace PFM {
	inline double derivativeX(const PFM::PeriodicDoublesLattice2D* field_ptr, const PFM::coordinate_t centerPoint) {
		
		return 0.5 *
				( field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y}) 
				  - field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y}) );
	}

	inline double derivativeY(const PFM::PeriodicDoublesLattice2D* field_ptr, const PFM::coordinate_t centerPoint) {

		return 0.5 *
				( field_ptr->getDataPoint({centerPoint.x, centerPoint.y + 1}) 
				  - field_ptr->getDataPoint({centerPoint.x, centerPoint.y - 1}) );
	}

	inline double laplacian5points(const PFM::PeriodicDoublesLattice2D* field_ptr, 
		                                      const PFM::coordinate_t centerPoint) {

	  return 0.25 * ( ( (-4) * field_ptr->getDataPoint({centerPoint.x, centerPoint.y}) )
				  + field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y})   
				  + field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y})
				  + field_ptr->getDataPoint({centerPoint.x, centerPoint.y + 1}) 
				  + field_ptr->getDataPoint({centerPoint.x, centerPoint.y - 1}) );
	}

	inline double laplacian5pointsAroundNeighCenter(const PFM::neighborhood9_t* neigh_ptr) {
	  return 0.25 * ( ( (-4) * neigh_ptr->data[4])
							  + neigh_ptr->data[1]   
							  + neigh_ptr->data[3]
							  + neigh_ptr->data[5] 
							  + neigh_ptr->data[7] 
					);
	}

	inline double laplacian9points(const PFM::PeriodicDoublesLattice2D* field_ptr, 
		                                      const PFM::coordinate_t centerPoint) {

	  return (1.0/3) * ( ( (-3) * field_ptr->getDataPoint({centerPoint.x, centerPoint.y}) )
					+ 0.5 * (field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y})   
								+ field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y})
								+ field_ptr->getDataPoint({centerPoint.x, centerPoint.y + 1}) 
								+ field_ptr->getDataPoint({centerPoint.x, centerPoint.y - 1}) 
							)
					+ 0.25 * (field_ptr->getDataPoint({centerPoint.x -1, centerPoint.y + 1})   
								+ field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y + 1})
								+ field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y - 1}) 
								+ field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y - 1}) 
							)
					);
					
	}

	inline double laplacian9pointsAroundNeighCenter(const PFM::neighborhood9_t* neigh_ptr) {
	  return (1.0/3) * ( ( (-3) * neigh_ptr->data[4])
								+ 0.5 * (neigh_ptr->data[1]  
										  + neigh_ptr->data[3]
										  + neigh_ptr->data[5] 
										  + neigh_ptr->data[7] 
								      )
								+ 0.25 * (neigh_ptr->data[0]  
										  + neigh_ptr->data[2]
										  + neigh_ptr->data[6] 
										  + neigh_ptr->data[8] 
								      )
					);
	}

	inline double trackAndReturnLaplacian9Contribution(const double weight, size_t index, 
		                                     PFM::neighborhood9_t* neigh_ptr, 
		                                     PFM::neighborhood9_t* neighContrib_ptr) {

		const double contrib = weight * neigh_ptr->data[index];
		neighContrib_ptr->data[index] += contrib;

		return contrib;
	}

	inline double laplacian9incrementingNeighborContributions(PFM::neighborhood9_t* neigh_ptr, 
							 				                  PFM::neighborhood9_t* neighContrib_ptr) {

		double lap = -3 * neigh_ptr->data[4];

		lap += trackAndReturnLaplacian9Contribution(0.5, 1, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.5, 3, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.5, 5, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.5, 7, neigh_ptr, neighContrib_ptr);

		lap += trackAndReturnLaplacian9Contribution(0.25, 0, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.25, 0, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.25, 6, neigh_ptr, neighContrib_ptr);
		lap += trackAndReturnLaplacian9Contribution(0.25, 8, neigh_ptr, neighContrib_ptr);

		return lap;
	}

	//Calculates the contribution from each neighbor in neigh and incrments it in the same index of neighContrib
	inline void calculateAndIncrementLaplacian9Contribution(PFM::neighborhood9_t* neigh_ptr, 
							 				                   PFM::neighborhood9_t* neighContrib_ptr) {

		neighContrib_ptr->data[1] += 0.5 * neigh_ptr->data[1];
		neighContrib_ptr->data[3] += 0.5 * neigh_ptr->data[3];
		neighContrib_ptr->data[5] += 0.5 * neigh_ptr->data[5];
		neighContrib_ptr->data[7] += 0.5 * neigh_ptr->data[7];

		neighContrib_ptr->data[0] += 0.25 * neigh_ptr->data[0];
		neighContrib_ptr->data[2] += 0.25 * neigh_ptr->data[2];
		neighContrib_ptr->data[6] += 0.25 * neigh_ptr->data[6];
		neighContrib_ptr->data[8] += 0.25 * neigh_ptr->data[8];

		return;
	}
}