#pragma once

//TODO: *in case* this is turned into a shared library, what should or shouldn't get a PFM_API appended?

//This exposes the data types for simulation parameters, configuration, checks and periodic 2D lattices.
//- simConfig_t deals with stuff like the size of the network, initial conditions, integration methods, etc;
//- simParameters_t holds the actual parameters, like dt, lambda, gamma, etc;
//- checkData_t holds runtime checks, like the avg density of the field or it's absolute change, etc;
//-- All of these have helper functions, specially to print out their data.
//- PeriodicDoublesLattice2D is a class which helps deal with and acess 2D periodic data (doubles).
//-- It also holds a public checkData_t and a private vector of them, for checks done on this lattice.
//--- Note that "checks done on this lattice" may or may not make sense depending on usage.

#include "PFM_API_enums.hpp"

#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/prng.hpp"

#include <assert.h>

#define ALL_CELLS_ID -1

namespace PFM {

	typedef struct parameterBounds_st {
		double minK = 0;
		double minA = 0;
		double minDt = 0;
		double maxK;
		double maxA;
		double maxDt;
	} parameterBounds_t;

	typedef struct simConfig_st {
		uint64_t stepsRan = 0;
		int cells = 0;
		int width = 0;
		int height = 0;
        double cellSeedValue = 0;
        uint64_t initialSeed = 0;
        simFuncEnum simulFunc = simFuncEnum::TOTAL_SIM_FUNCS;        
        initialConditions initialContidion = initialConditions::TOTAL_INITIAL_CONDS;
        double bias = -9999;
        integrationMethods method = integrationMethods::TOTAL_METHODS;
		bool perCellLayer = false;
		bool startPaused = false;

		std::chrono::system_clock::time_point epochTimeSimCall;
 
		inline uint32_t reducedSecondsSinceEpochOnSimCall() const {
			uint64_t seconds = 
				std::chrono::duration_cast<std::chrono::seconds>(epochTimeSimCall.time_since_epoch()).count();

			return (uint32_t)(seconds & UINT32_MAX);
		}
		
		inline std::string getSimDataString() const {
			std::string str = "";
			str += "Width: " + std::to_string(width) + "\n";
			str += "Height: " + std::to_string(height) + "\n";
			str += "Cells: " + std::to_string(cells) + "\n";
			str += "Bias: " + std::to_string(bias) + "\n";
			str += "PRNG Seed: " + std::to_string(initialSeed) + "\n";
			str += "Initial Cell Seed: " + std::to_string(cellSeedValue) + "\n";
			str += "Initial Condition: " + std::to_string((int)initialContidion) + "\n";
			str += "Method: " + std::to_string((int)method) + "\n";
			str += "Simulation: " + std::to_string((int)simulFunc) + "\n";
			str += "Steps: " + std::to_string(stepsRan) + "\n";
			str += "Per cell layer? " + std::to_string(perCellLayer) + "\n";
			str += "Start paused? " + std::to_string(startPaused);

			return str;
		}
	} simConfig_t;

	typedef struct simParameters_st {     
        double dt = -1;
		double lambda = -1;
		double gamma = -1;
        double A = -1;
		double k = -1;
		bool useAdaptativeDt = false;

		inline std::string getSimParamsString() const {
			std::string str = "";
			str += "dt: " + std::to_string(dt) + "\n";
			str += "Gamma (surface tension): " + std::to_string(gamma) + "\n";
			str += "Lambda (interface width): " + std::to_string(lambda) + "\n";
			str += "k: " + std::to_string(k) + "\n";
			str += "A: " + std::to_string(A) + "\n";
			str += "Adaptative dt? " + std::to_string(useAdaptativeDt);

			return str;
		}
	} simParameters_t;

	typedef struct coordinate_st {
		int x, y;

		bool operator==(coordinate_st const& otherCoord) const {
			return ((x == otherCoord.x) && (y == otherCoord.y));
		}
	} coordinate_t;

	//TODO: template for arbitrary sizes of (square) neighborhoods

	//in all methods, (x, y) is equivalent to (collum, row)
	typedef struct neighborhood9_st {
		double data[9];

		//TODO: add an indexOf method (how should it deal with bad input?)

		inline double getElement(int x, int y) {
			assert(x >= 0 && y >= 0);
			assert(y*3 + x < 9);
			return data[(y*3 + x)];
		}

		inline void setElement(int x, int y, double value) {
			assert(x >= 0 && y >= 0);
			assert(y*3 + x < 9);
			data[(y*3 + x)] = value;
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

		neighborhood9_st operator+=(neighborhood9_st const& otherNeigh) {
			for (size_t i = 0; i < 9; i++) { data[i] += otherNeigh.data[i]; }
			return *this;
		}

		neighborhood9_st operator-=(neighborhood9_st const& otherNeigh) {
			for (size_t i = 0; i < 9; i++) { data[i] -= otherNeigh.data[i]; }
			return *this;
		}

		neighborhood9_st operator+=(double const& scalar) {
			for (size_t i = 0; i < 9; i++) { data[i] += scalar; }
			return *this;
		}

		neighborhood9_st operator*=(double const& scalar) {
			for (size_t i = 0; i < 9; i++) { data[i] *= scalar; }
			return *this;
		}

	} neighborhood9_t;

	//in all methods, (x, y) is equivalent to (collum, row)
	typedef struct neighborhood25_st {
		double data[25];

		inline double getElement(int x, int y) {
			assert(x >= 0 && y >= 0);
			assert(y*5 + x < 25);
			return data[(y*5 + x)];
		}

		inline void setElement(int x, int y, double value) {
			assert(x >= 0 && y >= 0);
			assert(y*5 + x < 25);
			data[(y*5 + x)] = value;
		}

		inline double getCenter() const {
			return data[12];
		}

		inline void setCenter(double newValue) {
			data[12] = newValue;
		}

		inline void incrementCenter(double valueChange) {
			data[12] += valueChange;
		}

		neighborhood25_st operator+=(neighborhood25_st const& otherNeigh) {
			for (size_t i = 0; i < 25; i++) { data[i] += otherNeigh.data[i]; }
			return *this;
		}

		neighborhood25_st operator-=(neighborhood25_st const& otherNeigh) {
			for (size_t i = 0; i < 25; i++) { data[i] -= otherNeigh.data[i]; }
			return *this;
		}

		neighborhood25_st operator+=(double const& scalar) {
			for (size_t i = 0; i < 25; i++) { data[i] += scalar; }
			return *this;
		}

		neighborhood25_st operator*=(double const& scalar) {
			for (size_t i = 0; i < 25; i++) { data[i] *= scalar; }
			return *this;
		}

	} neighborhood25_t;

	typedef struct checkData_st {

		//TODO: rename stuff so it's clear what is reffering to the last step and what's about the last period
		//Also, what is per element or per unit time and what is not

		uint64_t step = 0;
		double absoluteChangeLastStep = 0;
		double sumOfsquaresOfChangesLastStep = 0;

		double totalChangeFromStart = 0;

		double lastDensity = 0;
		double densityChange = 0;
		double absoluteChange = 0;
		double sumOfsquaresOfChanges = 0;
		
		double lastDensityChange = 0;
		double lastAbsoluteChangePerElement = 0;
		double lastRmsAbsChange = 0;
		double lastAbsChangeStdDev = 0;

		double lastDt = 0;

		uint64_t stepsAtLastCheck = 0;
		int stepsDuringLastCheckPeriod = 0;
		double timeAtLastcheck = 0;
		double timeDuringLastCheckPeriod = 0;
		double totalTime = 0;
		double totalAbsoluteChangeSinceLastSave = 0;

		double remainingChangeSinceSaveOnLastCheck = 0;
		double referenceDt = 0;

		uint32_t substepsLastStep = 0;

		double totalAbsoluteChange = 0;
		
		simParameters_t parametersOnLastCheck;
		
		//TODO: also improve the naming of these
		inline void clearLastStepsChanges() { absoluteChangeLastStep = 0; sumOfsquaresOfChangesLastStep = 0; }
		inline void clearCurrentChanges() { densityChange = 0; absoluteChange = 0; sumOfsquaresOfChanges = 0; substepsLastStep = 0; }
		inline void zeroOut() { step = 0; totalChangeFromStart = 0; absoluteChangeLastStep = 0; 
							    sumOfsquaresOfChangesLastStep = 0; lastDensity = 0; densityChange = 0; 
								absoluteChange = 0; sumOfsquaresOfChanges = 0; timeAtLastcheck = 0; 
								timeDuringLastCheckPeriod = 0; lastDensityChange = 0; 
								lastAbsoluteChangePerElement = 0; lastRmsAbsChange = 0;
								lastAbsChangeStdDev = 0;  lastDt = 0; stepsAtLastCheck = 0;
		                        totalTime = 0; totalAbsoluteChangeSinceLastSave = 0;
		                        remainingChangeSinceSaveOnLastCheck = 0; referenceDt = 0; substepsLastStep = 0;
								totalAbsoluteChange = 0;
		}
		
		inline const std::string getChecksStr() const {
			
			assert(stepsDuringLastCheckPeriod > 0 && "Bad number of steps on last Check");
			
			const size_t precision = 10;
			char fmtValsBuffer[2*precision + 1];

			std::string str = "Checks @ step ";
			str += std::to_string(stepsAtLastCheck) + " (time: " + std::to_string(totalTime) 
				   + " steps during check: " + std::to_string(stepsDuringLastCheckPeriod) 
				   + " @avg dt: " + std::to_string(timeDuringLastCheckPeriod / stepsDuringLastCheckPeriod) + ")\n";
			str += "Density: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, lastDensity);
			str += fmtValsBuffer;
			str += "\n";
			str += "Change: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, lastDensityChange);
			str += fmtValsBuffer;
			str += ", absolute change / element per unit time: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, (lastAbsoluteChangePerElement / timeDuringLastCheckPeriod));
			str += fmtValsBuffer;
			str += " (since last save: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, totalAbsoluteChangeSinceLastSave);
			str += fmtValsBuffer;
			str += " - total: ";
			sprintf(fmtValsBuffer, "%.*f", 3, totalAbsoluteChange);
			str += fmtValsBuffer;
			str += ")\n";
			str += "Std Dev: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, lastAbsChangeStdDev / timeDuringLastCheckPeriod);
			str += fmtValsBuffer;
			str += " (RMS: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, lastRmsAbsChange / timeDuringLastCheckPeriod);
			str += fmtValsBuffer;
			str += ")\n";

			return str;
		}

		inline const std::string getParametersLastUpdateStr() const {

			const size_t precision = 10;
			char fmtValsBuffer[2*precision + 1];

			std::string str = "Parameters on the previous check:\nLambda: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, parametersOnLastCheck.lambda);
			str += fmtValsBuffer; 
			str += " | Gamma: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, parametersOnLastCheck.gamma);
			str += fmtValsBuffer; 
			str += " | dt: ";
			sprintf(fmtValsBuffer, "%.*f", (int)precision, parametersOnLastCheck.dt);
			str += fmtValsBuffer; 
			str += "\n";

			return str;
		}

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
		//Same, but 5x5.
		neighborhood25_t getNeighborhood25(coordinate_t centerPoint) const;
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void writeDataPoint(coordinate_t coordinate, double newValue);
		//WARNING: If an out-of-bounds element is asked, will ignore *with no warning*.
		void incrementDataPoint(coordinate_t coordinate, double changeInValue);
		
		//Adds to the checkData vector;
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

		const std::vector<checkData_t>* getCheckVectorConstPtr() const;

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

	//This holds 2 lattices to be used together, one for the current and another for the last value
	//If the checks vector is used, care should be taken to not confuse both (TODO: review this)
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