#pragma once

#include "PFM_data.hpp"

namespace PFM {
	inline double derivativeX(const PFM::PeriodicDoublesLattice2D* field_ptr, const PFM::coordinate_t centerPoint)
	{
		return 0.5 *
				( field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y}) 
				  - field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y}) );
	}

	inline double derivativeY(const PFM::PeriodicDoublesLattice2D* field_ptr, const PFM::coordinate_t centerPoint)
	{
		return 0.5 *
				( field_ptr->getDataPoint({centerPoint.x, centerPoint.y + 1}) 
				  - field_ptr->getDataPoint({centerPoint.x, centerPoint.y - 1}) );
	}

	inline double laplacian5points(const PFM::PeriodicDoublesLattice2D* field_ptr, 
		                                      const PFM::coordinate_t centerPoint)
	{
	  return 0.25 * ( ( (-4) * field_ptr->getDataPoint({centerPoint.x, centerPoint.y}) )
				  + field_ptr->getDataPoint({centerPoint.x + 1, centerPoint.y})   
				  + field_ptr->getDataPoint({centerPoint.x - 1, centerPoint.y})
				  + field_ptr->getDataPoint({centerPoint.x, centerPoint.y + 1}) 
				  + field_ptr->getDataPoint({centerPoint.x, centerPoint.y - 1}) );
	}
}