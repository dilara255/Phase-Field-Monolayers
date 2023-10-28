#pragma once

#include "miscStdHeaders.h"

//Multiple counters for different categories. 
//Can increment/decrement, get the counts (and the total count), or clear them.
//Will initialize the counters to zero on instantiation. Defaults to at least one category.
class MultipleCategoryCounter {

public:
	MultipleCategoryCounter(int totalCategories) {
		m_totalCategories = totalCategories;
		if (m_totalCategories < 1) {
			m_totalCategories = 1;
		}
		counter.reserve(m_totalCategories);
		for (int i = 0; i < m_totalCategories; i++) {
			counter.push_back(0);
		}
	}

	bool isCategoryIdValid(int categoryId) {
		if ((categoryId < 0) || (categoryId >= m_totalCategories)) {
			return false;
		}
		return true;
	}

	bool increment(int categoryId) {
		if(!isCategoryIdValid(categoryId)) {return false;}

		counter[categoryId]++;
		return true;
	}

	bool decrement(int categoryId) {
		if(!isCategoryIdValid(categoryId)) {return false;}

		counter[categoryId]--;
		return true;
	}

	int howMany(int categoryId) {
		if(!isCategoryIdValid(categoryId)) {return -1;}

		return counter[categoryId];
	}

	void printCounts() {
		printf("Count: 0: %d", counter[0]);
		for(int i = 1; i < m_totalCategories; i++){
			printf(", %d: %d", i, counter[i]);
		}
		puts("\n");
	}

	bool clearCounter(int categoryId) {
		if(!isCategoryIdValid(categoryId)) {return false;}
		
		counter[categoryId] = 0;
		return true;
	}

	void clearAll() {
		for(int i = 0; i < m_totalCategories; i++){
			counter[i] = 0;
		}		
	}

	int totalCount() {
		int total = 0;
		for(int i = 0; i < m_totalCategories; i++){
			total += counter[i];
		}	
		return total;
	}

private:
	int m_totalCategories;
	std::vector<int> counter;
};