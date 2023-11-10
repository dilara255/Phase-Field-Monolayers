#pragma once

#include "PFM_data.hpp"

inline double derivativeX(const PFM::PeriodicDoublesLattice2D& field, const PFM::coordinate_t centerPoint)
{
	return 0.5 *
			( field.getDataPoint({centerPoint.x + 1, centerPoint.y}) 
			  - field.getDataPoint({centerPoint.x - 1, centerPoint.y}) );
}

inline double derivativeY(const PFM::PeriodicDoublesLattice2D& field, const PFM::coordinate_t centerPoint)
{
	return 0.5 *
			( field.getDataPoint({centerPoint.x, centerPoint.y + 1}) 
			  - field.getDataPoint({centerPoint.x, centerPoint.y - 1}) );
}

inline double laplacian5points(const PFM::PeriodicDoublesLattice2D& field, const PFM::coordinate_t centerPoint)
{
  return ( (-4) * field.getDataPoint({centerPoint.x, centerPoint.y}) )
			+ field.getDataPoint({centerPoint.x + 1, centerPoint.y})   
			+ field.getDataPoint({centerPoint.x - 1, centerPoint.y})
			+ field.getDataPoint({centerPoint.x, centerPoint.y + 1}) 
			+ field.getDataPoint({centerPoint.x, centerPoint.y - 1});
}