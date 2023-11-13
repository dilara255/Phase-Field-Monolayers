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
}