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
				checks.density += initialData[i];
			}
			m_hasIntialized = true;
		}
		else {
			for (size_t i = 0; i < m_elements; i++) {
				m_data.push_back(0);
			}
		}
	}

	checks.density /= newDimensions.totalElements();
	
	//TODO: do I really want this printif here?
	printf("Field loaded with density %f\n", checks.density);

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
	if (index >= m_elements) { assert(false); return NAN; }
	else { return m_data[index]; }
}

//TODO: it might be helpful to optimize this a bit (but measures are still pending)
//TODO: Template for different neighborhood sizes?
neighborhood9_t PFM::PeriodicDoublesLattice2D::getNeighborhood(coordinate_t centerPoint) const {
	int centerX = centerPoint.x;
	int centerY = centerPoint.y;
	
	neighborhood9_t neigh = {
		getDataPoint({centerX - 1, centerY - 1}),
		getDataPoint({centerX + 0, centerY - 1}),
		getDataPoint({centerX + 1, centerY - 1}),

		getDataPoint({centerX - 1, centerY}),
		getDataPoint({centerX + 0, centerY}),
		getDataPoint({centerX + 1, centerY}),

		getDataPoint({centerX - 1, centerY + 1}),
		getDataPoint({centerX + 0, centerY + 1}),
		getDataPoint({centerX + 1, centerY + 1})
	};

	return neigh;
}

neighborhood25_t PFM::PeriodicDoublesLattice2D::getNeighborhood25(coordinate_t centerPoint) const {
	int centerX = centerPoint.x;
	int centerY = centerPoint.y;
	
	//SO UGLY :' )
	neighborhood25_t neigh = {
		
		getDataPoint({centerX - 2, centerY - 2}),
		getDataPoint({centerX - 1, centerY - 2}),
		getDataPoint({centerX + 0, centerY - 2}),
		getDataPoint({centerX + 1, centerY - 2}),
		getDataPoint({centerX + 2, centerY - 2}),

		getDataPoint({centerX - 2, centerY - 1}),
		getDataPoint({centerX - 1, centerY - 1}),
		getDataPoint({centerX + 0, centerY - 1}),
		getDataPoint({centerX + 1, centerY - 1}),
		getDataPoint({centerX + 2, centerY - 1}),

		getDataPoint({centerX - 2, centerY + 0}),
		getDataPoint({centerX - 1, centerY + 0}),
		getDataPoint({centerX + 0, centerY + 0}),
		getDataPoint({centerX + 1, centerY + 0}),
		getDataPoint({centerX + 2, centerY + 0}),

		getDataPoint({centerX - 2, centerY + 1}),
		getDataPoint({centerX - 1, centerY + 1}),
		getDataPoint({centerX + 0, centerY + 1}),
		getDataPoint({centerX + 1, centerY + 1}),
		getDataPoint({centerX + 2, centerY + 1}),

		getDataPoint({centerX - 2, centerY + 2}),
		getDataPoint({centerX - 1, centerY + 2}),
		getDataPoint({centerX + 0, centerY + 2}),
		getDataPoint({centerX + 1, centerY + 2}),
		getDataPoint({centerX + 2, centerY + 2})
	};

	return neigh;
}

void PFM::PeriodicDoublesLattice2D::writeDataPoint(coordinate_t coordinate, double newValue) {
	size_t index = indexFromPeriodicCoordinate(coordinate);
	if (index >= m_elements) { assert(false); return; }
	else { m_data[index] = newValue; return; }
}

void PFM::PeriodicDoublesLattice2D::incrementDataPoint(coordinate_t coordinate, double changeInValue) {
	size_t index = indexFromPeriodicCoordinate(coordinate);
	if (index >= m_elements) { assert(false); return; }
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

 bool PFM::PeriodicDoublesLattice2D::mirrorAllDataFrom(PeriodicDoublesLattice2D* otherField_ptr) {
	if(!m_hasIntialized || !m_hasAllocated) { return false; }
	if(!otherField_ptr->isInitialized() || !otherField_ptr->hasAllocated()) { return false; }
	if(m_data.size() != otherField_ptr->getNumberOfActualElements()) { return false; }
	
	for (size_t i = 0; i < m_data.size(); i++) {
		m_data[i] = otherField_ptr->getElement(i);
	}
	return true;
}

 const std::vector<checkData_t>* PFM::PeriodicDoublesLattice2D::getCheckVectorConstPtr() const {
	 return (const std::vector<checkData_t>*)&m_fieldChecks;
 }

 PFM::StringFmtBuffer::StringFmtBuffer(const size_t precision) {
	m_precision = precision;
	m_fmtValsBuffer_ptr = (char*)malloc(2 * precision + 1);
}

 PFM::StringFmtBuffer::~StringFmtBuffer() {
	free(m_fmtValsBuffer_ptr);
	m_fmtValsBuffer_ptr = nullptr;
}

template <typename T>
const char* PFM::StringFmtBuffer::getFormatedString(const T value) {
	write(value);
	return read();
}
template const char* PFM::StringFmtBuffer::getFormatedString<int64_t>(const int64_t param);
template const char* PFM::StringFmtBuffer::getFormatedString<uint64_t>(const uint64_t param);
template const char* PFM::StringFmtBuffer::getFormatedString<int>(const int param);
template const char* PFM::StringFmtBuffer::getFormatedString<unsigned int>(const unsigned int param);
template const char* PFM::StringFmtBuffer::getFormatedString<double>(const double param);
template const char* PFM::StringFmtBuffer::getFormatedString<float>(const float param);


const char* PFM::StringFmtBuffer::read() {
	if (m_fmtValsBuffer_ptr == nullptr) { return "Bad Buffer!\n"; }
	return (const char*)m_fmtValsBuffer_ptr;
}

//same as sprintif, except you don't pass the buffer
void PFM::StringFmtBuffer::write(const int64_t value) {
	if (m_fmtValsBuffer_ptr == nullptr) { return; }
	sprintf(m_fmtValsBuffer_ptr, "%*lld", (int)m_precision, value);
}

void PFM::StringFmtBuffer::write(const int value) {
	write(static_cast<int64_t>(value));
}

void PFM::StringFmtBuffer::write(const uint64_t value) {
	if (m_fmtValsBuffer_ptr == nullptr) { return; }
	sprintf(m_fmtValsBuffer_ptr, "%*llu", (int)m_precision, value);
}

void PFM::StringFmtBuffer::write(const unsigned int value) {
	write(static_cast<uint64_t>(value));
}

void PFM::StringFmtBuffer::write(const double value) {
	if (m_fmtValsBuffer_ptr == nullptr) { return; }
	sprintf(m_fmtValsBuffer_ptr, "%.*f", (int)m_precision, value);
}

void PFM::StringFmtBuffer::write(const float value) {
	write(static_cast<double>(value));
}
 