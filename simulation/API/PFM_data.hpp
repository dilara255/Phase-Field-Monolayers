#pragma once

//TODO: *in case* this is turned into a shared library, what should or shouldn't get a PFM_API appended?

#include <vector>

namespace PFM {

	typedef struct coordinate_str {
		int x, y;
	} coordinate_t;

	typedef struct fieldDimensions_str {
		size_t width, height;
		size_t totalElements() const;
	} fieldDimensions_t;

	class PeriodicDoublesLattice2D {
	public:
		//If invalid dimensions are passed, no memory is allocated. 
		//If no initialData is passed, no data is initialized, but memory is allocated.
		//If the size of the data passed doesn't fit the dimensions *exactly*, it's ignored.
		PeriodicDoublesLattice2D(fieldDimensions_t newDimensions, std::vector<double> initialData = {});

		bool hasAllocated() const;
		bool isInitialized() const;		
		size_t getNumberOfAllocatedElements() const;
		size_t getNumberOfActualElements() const;
		fieldDimensions_t getFieldDimensions() const;

		//WARNING: *WON'T* check wether the field is initialized: will return junk if it's not!
		//If an out-of-bounds element is asked, will return NAN.
		double getDataPoint(coordinate_t coordinate) const;
		double getElement(size_t index) const;
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void writeDataPoint(coordinate_t coordinate, double newValue);
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void incrementDataPoint(coordinate_t coordinate, double changeInValue);
		
		//Allocates a new buffer and passes the data. *dimensions_ptr will hold the dimensions of the new buffer.
		//The new buffer becomes a responsability of the caller.
		//Will return null if not yet allocated or initialized or if the new allocation fails.
		double* createCopyOfAllData(fieldDimensions_t* dimensions_ptr) const;
		
		//If elements in newData is different from the (actual) size of the field, returns false.
		//Otherwise, accepts the new data and marks the field as initialized.
		bool rewriteAllData(std::vector<double> newData);

		//If elements in newData > the size of the field after the starting point, returns false.
		//If the field is not yet initialized or allocated, also returns false (and doesn't accept new data).
		bool acceptBulkData(coordinate_t startingPoint, std::vector<double> newData);
	
	private:
		std::vector<double> m_data;
		bool m_hasAllocated;
		bool m_hasIntialized;
		fieldDimensions_t m_dimensions;
		size_t m_elements;

		size_t indexFromBoundedCoordinate(coordinate_t coordinate) const;
		size_t indexFromPeriodicCoordinate(coordinate_t coordinate) const;
	};
}