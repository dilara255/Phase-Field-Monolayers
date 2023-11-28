#pragma once

//TODO: *in case* this is turned into a shared library, what should or shouldn't get a PFM_API appended?

#include "PFM_API_enums.hpp"

#include <vector>
#include <assert.h>
#include <memory>

#define ALL_CELLS_ID -1

namespace PFM {

	typedef struct simData_st {
		int stepsRan = 0;
		int cells = 0;
        double lastCellSeedValue = 0;
        uint64_t initialSeed = 0;
        simFuncEnum lastSimulFuncUsed = simFuncEnum::TOTAL_SIM_FUNCS;        
        double lastK = -1;
        double lastA = -1;
        double lastDT = -1;
        initialConditions lastInitialContidion = initialConditions::TOTAL_INITIAL_CONDS;
        double lastBias = -9999;
        integrationMethods lastMethod = integrationMethods::TOTAL_METHODS;
	} simData_t;

	typedef struct coordinate_st {
		int x, y;
	} coordinate_t;

	typedef struct neighborhood9_st {
		double data[9];

		inline double getElement(int x, int y) {
			assert(x >= 0 && y >= 0);
			assert(y*3 + x < 9);
			return data[(y*3 + x)];
		}

		inline double getCenter() const {
			return data[4];
		}

		inline void setCenter(double newValue) {
			data[4] = newValue;
		}

		inline void incrementCenter(double valueChange) {
			data[4] += valueChange;
		}

	} neighborhood9_t;

	typedef struct checkData_st {
		int step = 0;
		double lastDensity = 0;
		double densityChange = 0;
		double absoluteChange = 0;

		inline void clearChanges() {densityChange = 0; absoluteChange = 0;}
		inline void zeroOut() { step = 0; lastDensity = 0; densityChange = 0; absoluteChange = 0; }
	} checkData_t;

	typedef struct fieldDimensions_st {
		size_t width, height;
		size_t totalElements() const;
	} fieldDimensions_t;

	class PeriodicDoublesLattice2D {
	public:
		//If invalid dimensions are passed, no memory is allocated. 
		//If no initialData is passed, no data is initialized, but memory is allocated.
		//If the size of the data passed doesn't fit the dimensions *exactly*, the field is filled with 0s.
		PeriodicDoublesLattice2D(fieldDimensions_t newDimensions, int cellID = ALL_CELLS_ID, 
			                                           std::vector<double> initialData = {});

		bool hasAllocated() const;
		bool isInitialized() const;		
		size_t getNumberOfAllocatedElements() const;
		size_t getNumberOfActualElements() const;
		fieldDimensions_t getFieldDimensions() const;

		//WARNING: *WON'T* check wether the field is initialized: will return junk if it's not!
		//If an out-of-bounds element is asked, will return NAN.
		double getDataPoint(coordinate_t coordinate) const;
		double getElement(size_t index) const;
		//Gets a structure with the 3x3 neighboorhod of the centerpoint. 
		neighborhood9_t getNeighborhood(coordinate_t centerPoint) const;
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void writeDataPoint(coordinate_t coordinate, double newValue);
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void incrementDataPoint(coordinate_t coordinate, double changeInValue);
		
		//Adds a density and it's step to a vector with the previous field densities;
		void addFieldCheckData(checkData_t checkData);

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

		//If the fields don't have the exact same size, returns false.
		//If either of the field are not yet initialized or allocated, also returns false.
		bool mirrorAllDataFrom(PeriodicDoublesLattice2D* otherField_ptr);

		//To hold check data for the field. Is initialized on field initialization
		//Other than that, is of responsability of the user of the class
		checkData_t checks;

	private:
		std::vector<double> m_data;

		bool m_hasAllocated;
		bool m_hasIntialized;
		fieldDimensions_t m_dimensions;
		size_t m_elements;

		std::vector<checkData_t> m_fieldChecks;
		
		size_t indexFromBoundedCoordinate(coordinate_t coordinate) const;
		size_t indexFromPeriodicCoordinate(coordinate_t coordinate) const;
	};

	//This holds 2 lattices to be used in conjunction, one for the current and another for the last value
	class CurrentAndLastPerioricDoublesLattice2D {
	public:
		CurrentAndLastPerioricDoublesLattice2D(fieldDimensions_t newDimensions, int cellID = ALL_CELLS_ID, 
			                                                            std::vector<double> initialData = {});
		
		PeriodicDoublesLattice2D* getPointerToCurrent();
		PeriodicDoublesLattice2D* getPointerToLast();

		void rotatePointers(); 

		void releaseFields();

		//TODO: write method to write into both fields

	private:
		std::unique_ptr<PeriodicDoublesLattice2D> m_lattice1_ptr;
		std::unique_ptr<PeriodicDoublesLattice2D> m_latticeA_ptr;
		
		bool m_isLattice1current = true;
	};
}