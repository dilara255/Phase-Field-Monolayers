#include "PFM_data.hpp"

#include <cassert>

using namespace PFM;

size_t PFM::fieldDimensions_st::totalElements() const {
	return width * height;
}

PFM::CurrentAndLastPerioricDoublesLattice2D::CurrentAndLastPerioricDoublesLattice2D(fieldDimensions_t newDimensions, 
																				int cellID,
																				std::vector<double> initialData) {

	m_latticeA_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
		new PeriodicDoublesLattice2D(newDimensions, cellID, initialData)
	);
	m_lattice1_ptr = std::unique_ptr<PeriodicDoublesLattice2D>(
		new PeriodicDoublesLattice2D(newDimensions, cellID, initialData)
	);
}
		
PFM::PeriodicDoublesLattice2D* PFM::CurrentAndLastPerioricDoublesLattice2D::getPointerToCurrent() {
	if(m_isLattice1current) { return m_lattice1_ptr.get(); }
	else { return m_latticeA_ptr.get(); }
}

PFM::PeriodicDoublesLattice2D* PFM::CurrentAndLastPerioricDoublesLattice2D::getPointerToLast() {
	if(m_isLattice1current) { return m_latticeA_ptr.get(); }
	else { return m_lattice1_ptr.get(); }
}

void PFM::CurrentAndLastPerioricDoublesLattice2D::rotatePointers() {
	m_isLattice1current = !m_isLattice1current;
}

void PFM::CurrentAndLastPerioricDoublesLattice2D::releaseFields() {
	if(m_latticeA_ptr != NULL) { m_latticeA_ptr.release(); } 
	if(m_lattice1_ptr != NULL) { m_lattice1_ptr.release(); } 
}

PFM::PeriodicDoublesLattice2D::PeriodicDoublesLattice2D(fieldDimensions_t newDimensions, int cellID, 
	                                                                std::vector<double> initialData) {
	m_hasAllocated = false;
	m_hasIntialized = false;

	if (newDimensions.width <= 0 || newDimensions.height <= 0) { 
		m_elements = 0;
		//but just to avoid divisions by zero:
		m_dimensions.width = 1;
		m_dimensions.height = 1;
		
		return; 
	}

	m_dimensions = newDimensions;
	m_elements = m_dimensions.totalElements();
	m_data.reserve(m_elements);
	m_hasAllocated = true;

	if(cellID == ALL_CELLS_ID) {
		if (initialData.size() == m_elements) {
			for (size_t i = 0; i < m_elements; i++) {
				m_data.push_back(initialData[i]);
			}
			m_hasIntialized = true;
		}
	}

	return;
}

bool PFM::PeriodicDoublesLattice2D::hasAllocated() const {
	return m_hasAllocated;
}

bool PFM::PeriodicDoublesLattice2D::isInitialized() const {
	return m_hasIntialized;
}

 size_t PFM::PeriodicDoublesLattice2D::getNumberOfAllocatedElements() const {
	return m_data.capacity();
}

 size_t PFM::PeriodicDoublesLattice2D::getNumberOfActualElements() const {
	return m_data.size();
}

fieldDimensions_t PFM::PeriodicDoublesLattice2D::getFieldDimensions() const {
	return m_dimensions;
}

size_t PFM::PeriodicDoublesLattice2D::indexFromBoundedCoordinate(coordinate_t coordinate) const {
	assert(coordinate.x >= 0 && coordinate.y >= 0 && 
		   coordinate.x < (int)m_dimensions.width && coordinate.y < (int)m_dimensions.height);

	return coordinate.x + (coordinate.y * m_dimensions.width);
}

size_t PFM::PeriodicDoublesLattice2D::indexFromPeriodicCoordinate(coordinate_t coordinate) const {
	//% with negative is implementation defined, so, to avoid ub we first take the values as positve:
	size_t effectiveX = std::abs(coordinate.x) % m_dimensions.width;
	size_t effectiveY = std::abs(coordinate.y) % m_dimensions.height;

	//and then adjust in case they were negative (the "2*" is meant to "clean" the positive and then "discount"):
	effectiveY += (coordinate.y < 0) * (m_dimensions.height - (2*effectiveY) );
	effectiveX += (coordinate.x < 0) * (m_dimensions.width - (2*effectiveX) );

	assert(effectiveX >= 0 && effectiveY >= 0 && 
		   effectiveX < m_dimensions.width && effectiveY < m_dimensions.height);

	return effectiveX + (effectiveY * m_dimensions.width);
}

double PFM::PeriodicDoublesLattice2D::getDataPoint(coordinate_t coordinate) const {
	size_t index = indexFromPeriodicCoordinate(coordinate);
	return getElement(index);
}

double PFM::PeriodicDoublesLattice2D::getElement(size_t index) const {
	if (index >= m_elements) { return NAN; }
	else { return m_data[index]; }
}

//TODO: it might be helpful to optimize this a bit (but measures are still pending)
neighborhood9_t PFM::PeriodicDoublesLattice2D::getNeighborhood(coordinate_t centerPoint) const {
	int centerX = centerPoint.x;
	int centerY = centerPoint.y;
	
	neighborhood9_t neigh = {
		getDataPoint({centerX - 1, centerY - 1}),
		getDataPoint({centerX, centerY - 1}),
		getDataPoint({centerX + 1, centerY - 1}),

		getDataPoint({centerX - 1, centerY}),
		getDataPoint({centerX, centerY}),
		getDataPoint({centerX + 1, centerY}),

		getDataPoint({centerX - 1, centerY + 1}),
		getDataPoint({centerX, centerY + 1}),
		getDataPoint({centerX + 1, centerY + 1})
	};

	return neigh;
}

void PFM::PeriodicDoublesLattice2D::writeDataPoint(coordinate_t coordinate, double newValue) {
	size_t index = indexFromPeriodicCoordinate(coordinate);
	if (index >= m_elements) { return; }
	else { m_data[index] = newValue; return; }
}

void PFM::PeriodicDoublesLattice2D::incrementDataPoint(coordinate_t coordinate, double changeInValue) {
	size_t index = indexFromPeriodicCoordinate(coordinate);
	if (index >= m_elements) { return; }
	else { m_data[index] += changeInValue; return; }
}

void PFM::PeriodicDoublesLattice2D::addFieldCheckData(PFM::checkData_t checkData) {
	m_fieldChecks.push_back(checkData);
}

 double* PFM::PeriodicDoublesLattice2D::createCopyOfAllData(fieldDimensions_t* dimensions_ptr) const {
	*dimensions_ptr = m_dimensions;
	if (!m_hasAllocated || !m_hasIntialized) { return NULL; }

	double* data_ptr = (double*)malloc( (sizeof(double)) * m_elements);
	if(data_ptr == NULL) { return NULL; }

	for (size_t i = 0; i < m_elements; i++) {
		data_ptr[i] = m_data[i];
	}
	return data_ptr;
}

 bool PFM::PeriodicDoublesLattice2D::rewriteAllData(std::vector<double> newData) {
	if(newData.size() != m_data.size()) { return false; }

	for (size_t i = 0; i < m_elements; i++) {
		m_data[i] = newData[i];
	}

	return true;
}

 bool PFM::PeriodicDoublesLattice2D::acceptBulkData(coordinate_t startingPoint, 
	                                                     std::vector<double> newData) {
	if(!m_hasIntialized || !m_hasAllocated) { return false; }

	size_t startingIndex = indexFromPeriodicCoordinate(startingPoint);
	if(startingIndex + newData.size() > m_elements) { return false; }
	
	for (size_t i = 0; i < newData.size(); i++) {
		m_data[startingIndex + i] = newData[i];
	}
	return true;
}