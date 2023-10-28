#pragma once

#include "miscStdHeaders.h"

#include "prng.hpp"
#include "timeHelpers.hpp"
#include "flagFields.hpp"

//TODO: move declarations to AS_testsAPIT and move the definitions to cpp file
namespace AZ{

	static bool specificFlagField32Test(FlagField32* flags, uint32_t testField,
  		     int testBit, int initialExpectedOn, int testBitInitialExpectation, 
		                                                      bool log = false) {
		int testsDone = 0;

		flags->loadField(testField);
		bool result = (flags->isBitOn(testBit) == (bool)testBitInitialExpectation);
		int expectedOn = initialExpectedOn;
		result &= (flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
					flags->howManyAreOn(), expectedOn, flags->getFlagsAsString().c_str());
		}

		flags->setBitOn(testBit);
		result &= (flags->isBitOn(testBit) == true);
		expectedOn = (initialExpectedOn + (1 - testBitInitialExpectation));
		result &= 
			(flags->howManyAreOn() == expectedOn);

		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
					flags->howManyAreOn(), expectedOn, flags->getFlagsAsString().c_str());
		}

		flags->setBitOff(testBit);
		result &= (flags->isBitOn(testBit) == false);
		expectedOn = (initialExpectedOn - testBitInitialExpectation);
		result &= (flags->howManyAreOn() == expectedOn);

		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
					flags->howManyAreOn(), expectedOn, flags->getFlagsAsString().c_str());
		}

		flags->setBitOff(testBit);
		result &= (flags->isBitOn(testBit) == false);
		expectedOn = (initialExpectedOn - testBitInitialExpectation);
		result &= (flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
					flags->howManyAreOn(), expectedOn, flags->getFlagsAsString().c_str());
		}

		flags->setBitOn(testBit);
		result &= (flags->isBitOn(testBit) == true);
		expectedOn = (initialExpectedOn + (1 - testBitInitialExpectation));
		result &= 
			(flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
					flags->howManyAreOn(), expectedOn, flags->getFlagsAsString().c_str());
		}

		return result;
	}

	static bool specificFlagField128Test(FlagField128* flags, uint32_t testField, int blockTo,
						    int testBit, int initialExpectedOn, int testBitInitialExpectation,
   		                                                                     bool log = false) {
		int testsDone = 0;

		flags->loadField(testField, blockTo);
		bool result = (flags->isBitOn(testBit) == (bool)testBitInitialExpectation);
		int expectedOn = initialExpectedOn;
		result &= (flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;			
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
				flags->howManyAreOn(), expectedOn, flags->getFlagsAsString(blockTo).c_str());
		}

		flags->setBitOn(testBit);
		result &= (flags->isBitOn(testBit) == true);
		expectedOn = (initialExpectedOn + (1 - testBitInitialExpectation));
		result &= 
			(flags->howManyAreOn() == expectedOn);

		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
				flags->howManyAreOn(), expectedOn, flags->getFlagsAsString(blockTo).c_str());
		}

		flags->setBitOff(testBit);
		result &= (flags->isBitOn(testBit) == false);
		expectedOn = (initialExpectedOn - testBitInitialExpectation);
		result &= (flags->howManyAreOn() == expectedOn);

		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
				flags->howManyAreOn(), expectedOn, flags->getFlagsAsString(blockTo).c_str());
		}

		flags->setBitOff(testBit);
		result &= (flags->isBitOn(testBit) == false);
		expectedOn = (initialExpectedOn - testBitInitialExpectation);
		result &= (flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
				flags->howManyAreOn(), expectedOn, flags->getFlagsAsString(blockTo).c_str());
		}

		flags->setBitOn(testBit);
		result &= (flags->isBitOn(testBit) == true);
		expectedOn = (initialExpectedOn + (1 - testBitInitialExpectation));
		result &= 
			(flags->howManyAreOn() == expectedOn);
		
		if(log) {
			testsDone++;
			printf("\n %d: %d (on: %d, expected: %d, flags: %s)", testsDone, (int)result, 
				flags->howManyAreOn(), expectedOn, flags->getFlagsAsString(blockTo).c_str());
		}

		return result;
	}

	//Creates one flagField of each kind, checks:
	//-Changes and retrieval (all modes for each bank); and
	//-Count of on flags, along the way.
	//Returns false if any test fails.
	static bool testFlagFields(bool log = false) {
		
		LOG_DEBUG("Will test Flag Field functionality\n", 1); 
		GETCHAR_PAUSE;

		uint32_t full = 0b11111111111111111111111111111111;
		uint32_t empty = 0b00000000000000000000000000000000;
		uint32_t some = 0b10101010101010101010101010101010;
		uint32_t other = 0b01010101010101010101010101010101;

		int expectedOnFull = 32;
		int expectedOnEmpty = 0;
		int expectedOnSome = 16;
		int expectedOnOther = 16;

		int full_0 = 1;
		int empty_0 = 0;
		int some_0 = 0;
		int other_0 = 1;

		FlagField32 testField32;

		int bitToTest = 0;

		bool result = specificFlagField32Test(&testField32, full, bitToTest, expectedOnFull, 
			                                                                   full_0, log);
		if(log) { puts("\n----------------"); }

		int expectOn = testField32.howManyAreOn();
		uint32_t field = testField32.getField();
		bitToTest = 10;
		int testBitExpected = testField32.isBitOn(bitToTest);
		result &= specificFlagField32Test(&testField32, field, bitToTest, 
			                                   expectOn, testBitExpected, log);
		if(log) { puts("\n----------------"); }

		expectOn = testField32.howManyAreOn();
		field = testField32.getField();
		bitToTest = 31;
		testBitExpected = testField32.isBitOn(bitToTest);
		result &= specificFlagField32Test(&testField32, field, bitToTest, 
			                                   expectOn, testBitExpected, log);
		if(log) { puts("\n----------------"); }

		bitToTest = 0;
		result &= specificFlagField32Test(&testField32, empty, bitToTest, expectedOnEmpty, 
			                                                                 empty_0, log);
		if(log) { puts("\n----------------"); }

		expectOn = testField32.howManyAreOn();
		field = testField32.getField();
		bitToTest = 10;
		testBitExpected = testField32.isBitOn(bitToTest);
		result &= specificFlagField32Test(&testField32, field, bitToTest, 
			                                   expectOn, testBitExpected, log);
		if(log) { puts("\n----------------"); }

		expectOn = testField32.howManyAreOn();
		field = testField32.getField();
		bitToTest = 31;
		testBitExpected = testField32.isBitOn(bitToTest);
		result &= specificFlagField32Test(&testField32, field, bitToTest, 
			                                   expectOn, testBitExpected, log);
		if(log) { puts("\n----------------"); }
		if(log) { puts("\n----------------"); }

		FlagField128 testfield128;

		int blockToTest = 0;
		bitToTest = 0;
		result &= specificFlagField128Test(&testfield128, full, blockToTest, bitToTest,
													       expectedOnFull, full_0, log);
		if(log) { puts("\n----------------"); }

		field = empty;
		expectOn = testfield128.howManyAreOn() + expectedOnEmpty;
		blockToTest = 1;
		bitToTest = 0 + 32*blockToTest;
		result &= specificFlagField128Test(&testfield128, field, blockToTest, bitToTest,
																 expectOn, empty_0, log);
		if(log) { puts("\n----------------"); }

		field = some;
		expectOn = testfield128.howManyAreOn() + expectedOnSome;
		blockToTest = 2;
		bitToTest = 0 + 32*blockToTest;
		result &= specificFlagField128Test(&testfield128, field, blockToTest, bitToTest,
																  expectOn, some_0, log);
		if(log) { puts("\n----------------"); }

		field = other;
		expectOn = testfield128.howManyAreOn() + expectedOnOther;
		blockToTest = 3;
		bitToTest = 0 + 32*blockToTest;
		result &= specificFlagField128Test(&testfield128, field, blockToTest, bitToTest,
																 expectOn, other_0, log);
		if(log) { puts("\n----------------"); }

		field = testfield128.getField(3);
		expectOn = testfield128.howManyAreOn();
		blockToTest = 3;
		bitToTest = 31 + 32*blockToTest;
		testBitExpected = testfield128.isBitOn(bitToTest);
		result &= specificFlagField128Test(&testfield128, field, blockToTest, bitToTest,
														 expectOn, testBitExpected, log);
		
		if(!result) { LOG_ERROR("Something went wrong!"); }
		else if(log) { LOG_INFO("Passed"); }
		
		return result;
	}

	//Draws howManyToDraw prns into a vector. Tests:
    //Average, amount of zeroes, periods 2, 4, 8 and 16.
    //Returns true if all tests pass. Defaults to 512k draws, minimum 16.
	//If log == true, logs to output (eg, console) the time per draw
	static bool testDraw4spcg32s(bool log = true, uint32_t howManyTuplesToDraw = 128000) {
		
		LOG_DEBUG("Will test PRN system: draw 4 at a time\n", 1); 
		GETCHAR_PAUSE;

		uint64_t totalToDraw = (uint64_t)howManyTuplesToDraw * DRAW_WIDTH;
		if (totalToDraw >= (std::numeric_limits<uint32_t>::max())) {
			LOG_ERROR("Too many tuples were asked - might overflow: aborting test");
			return false;
		}

		const int mininumTotalDrawn = 16;

		if (howManyTuplesToDraw < (mininumTotalDrawn / DRAW_WIDTH)) {
			if (log) {
				puts("Too few tuples for test to run properly, aborting with false");
			}
			return false;
		}

		//Draws howManyToDraw prns into a vector, timing it:
		std::vector<uint32_t> drawn;
		uint32_t totalDrawn = DRAW_WIDTH * howManyTuplesToDraw;
		drawn.reserve(totalDrawn);
		for (uint32_t i = 0; i < totalDrawn; i++) {
			drawn.push_back(0);
		}

		uint64_t seeds[DRAW_WIDTH];
		seeds[0] = DEFAULT_PRNG_SEED0;
		seeds[1] = DEFAULT_PRNG_SEED1;
		seeds[2] = DEFAULT_PRNG_SEED2;
		seeds[3] = DEFAULT_PRNG_SEED3;

		auto startTime = std::chrono::steady_clock::now();	
		
		for (uint32_t i = 0; i < howManyTuplesToDraw; i++) {
			
			draw4spcg32s(&seeds[0],&seeds[1],&seeds[2],&seeds[3],
							&drawn[i*DRAW_WIDTH], &drawn[(i*DRAW_WIDTH) + 1],
							&drawn[(i*DRAW_WIDTH) + 2], &drawn[(i*DRAW_WIDTH) + 3]);
		}

		auto endTime = std::chrono::steady_clock::now();
		auto timeElapsed = endTime - startTime;
		double nanosPerDraw = (double)timeElapsed.count()/totalDrawn;

		//Tests amount of zeroes:

		int zeroes = 0;
		for (uint32_t i = 0; i < totalDrawn; i++) {	
			zeroes += (drawn[i] == 0);
		}

		double expectedZeroes = (double)totalDrawn/UINT32_MAX;
		double proportionalDifference = fabs((zeroes - expectedZeroes))/expectedZeroes;
		int64_t graceZeroes = 1;
		double maxProportionalDifference = graceZeroes/expectedZeroes;

		bool result = (proportionalDifference <= maxProportionalDifference);
		
		//Tests difference in periods 2, 4, 8 and 16.
		//For values in (0,1), the expected difference is given by:
		//E[|X - Y|] = P(Y > X)E[Y - X | Y > X] + P(Y <= X)E[X - Y | X > Y] ,
		//which turns out to be 1/3. The expected variance is 1/18.
		double expectedAverageAbsoluteDifference = 1.0/3;

		//Will take to (0,1) first
		std::vector<double> drawnNormalized;
		drawnNormalized.reserve(totalDrawn);

		double inverseOfUint32max = 1.0/UINT32_MAX;
		for (uint32_t i = 0; i < totalDrawn; i++) {	
			drawnNormalized.push_back(drawn[i]*inverseOfUint32max);
		}

		double inverseOfTotalNumbersDrawn = 1.0/(totalDrawn - mininumTotalDrawn);
		double average = 0;
		double averageAbsoluteDifferencePeriod2 = 0;
		double averageAbsoluteDifferencePeriod4 = 0;
		double averageAbsoluteDifferencePeriod8 = 0;
		double averageAbsoluteDifferencePeriod16 = 0;
		for (uint32_t i = mininumTotalDrawn; i < totalDrawn; i++) {	
			average += drawnNormalized[i]*inverseOfTotalNumbersDrawn;
			averageAbsoluteDifferencePeriod2 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-2])*inverseOfTotalNumbersDrawn;
			averageAbsoluteDifferencePeriod4 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-4])*inverseOfTotalNumbersDrawn;
			averageAbsoluteDifferencePeriod8 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-8])*inverseOfTotalNumbersDrawn;
			averageAbsoluteDifferencePeriod16 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-16])*inverseOfTotalNumbersDrawn;
		}

		double margin = sqrt(inverseOfTotalNumbersDrawn);
		double differenceFromExpectedPeriod2 = 
			fabs(averageAbsoluteDifferencePeriod2 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod4 = 
			fabs(averageAbsoluteDifferencePeriod4 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod8 = 
			fabs(averageAbsoluteDifferencePeriod8 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod16 = 
			fabs(averageAbsoluteDifferencePeriod16 - expectedAverageAbsoluteDifference);

		result &= (differenceFromExpectedPeriod2 <= margin);
		result &= (differenceFromExpectedPeriod4 <= margin);
		result &= (differenceFromExpectedPeriod8 <= margin);
		result &= (differenceFromExpectedPeriod16 <= margin);

		//Log and return result:

		if (log) {
			printf("\nPRNG test (draw4spcg32s): %d drawn (normalized avg: %f)\n",
				                                               totalDrawn, average);
			printf("\tTiming: %f nanos per prn\n", nanosPerDraw);
			printf("\tZeroes: %d (expected: %f)\n", zeroes, expectedZeroes);
			printf("\tAvg. differences expected: %f (margin used: %f). Found:\n",
								expectedAverageAbsoluteDifference, margin);
			printf("\t\tPeriod 2: %f; Period 4: %f, Period 8: %f, Period 16: %f\n",
				averageAbsoluteDifferencePeriod2, averageAbsoluteDifferencePeriod4,
			   averageAbsoluteDifferencePeriod8, averageAbsoluteDifferencePeriod16);
		}

		if(!result) { LOG_ERROR("Something went wrong!"); }
		else if(log) { LOG_INFO("Passed"); }

		return result;
	}

    //Draws howManyToDraw prns into a vector. Tests:
    //Average, amount of zeroes, periods 2, 4, 8 and 16
    //Returns true if all tests pass. Defaults to 128k draws
	//If log == true, logs to output (eg, console) the time per draw
	static bool testDraw1spcg32(bool log = true, uint32_t howManyToDraw = 128000) {
			const int mininumTotalDrawn = 16;

		LOG_DEBUG("Will test PRN system: draw 1 at a time\n", 1); 
		GETCHAR_PAUSE;

		if (howManyToDraw < mininumTotalDrawn) {
			if (log) {
				LOG_ERROR("Too few tuples for test to run properly, aborting with false");
			}
			return false;
		}

		//Draws howManyToDraw prns into a vector, timing it:
		std::vector<uint32_t> drawn;
		uint32_t totalDrawn = howManyToDraw;
		drawn.reserve(totalDrawn);
		for (uint32_t i = 0; i < totalDrawn; i++) {	
			drawn.push_back(0);
		}

		uint64_t seed = DEFAULT_PRNG_SEED0;

		auto startTime = std::chrono::steady_clock::now();

		for (uint32_t i = 0; i < totalDrawn; i++) {	
			drawn[i] = draw1spcg32(&seed);
		}

		auto endTime = std::chrono::steady_clock::now();
		auto timeElapsed = endTime - startTime;
		double nanosPerDraw = (double)timeElapsed.count()/totalDrawn;

		//Tests amount of zeroes:
		int zeroes = 0;
		for (uint32_t i = 0; i < totalDrawn; i++) {	
			zeroes += (drawn[i] == 0);
		}

		double expectedZeroes = (double)totalDrawn/UINT32_MAX;
		double proportionalDifference = fabs((zeroes - expectedZeroes))/expectedZeroes;
		int64_t graceZeroes = 1;
		double maxProportionalDifference = graceZeroes/expectedZeroes;

		bool result = (proportionalDifference <= maxProportionalDifference);
		
		//Tests difference in periods 2, 4, 8 and 16.
		//For values in (0,1), the expected difference is given by:
		//E[|X - Y|] = P(Y > X)E[Y - X | Y > X] + P(Y <= X)E[X - Y | X > Y] ,
		//which turns out to be 1/3. The expected variance is 1/18.
		double expectedAverageAbsoluteDifference = 1.0/3;

		//Will take to (0,1) first
		std::vector<double> drawnNormalized;
		drawnNormalized.reserve(totalDrawn);

		double inverseOfUint32max = 1.0/UINT32_MAX;
		for (uint32_t i = 0; i < totalDrawn; i++) {	
			drawnNormalized.push_back(drawn[i]*inverseOfUint32max);
		}

		double inverseOfTotalNumbersTested = 1.0/(totalDrawn - mininumTotalDrawn);
		double avarege = 0;
		double averageAbsoluteDifferencePeriod2 = 0;
		double averageAbsoluteDifferencePeriod4 = 0;
		double averageAbsoluteDifferencePeriod8 = 0;
		double averageAbsoluteDifferencePeriod16 = 0;

		for (uint32_t i = mininumTotalDrawn; i < totalDrawn; i++) {	
			avarege += drawnNormalized[i]*inverseOfTotalNumbersTested;
			averageAbsoluteDifferencePeriod2 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-2])*inverseOfTotalNumbersTested;
			averageAbsoluteDifferencePeriod4 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-4])*inverseOfTotalNumbersTested;
			averageAbsoluteDifferencePeriod8 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-8])*inverseOfTotalNumbersTested;
			averageAbsoluteDifferencePeriod16 +=
				fabs(drawnNormalized[i] - drawnNormalized[i-16])*inverseOfTotalNumbersTested;
		}

		double margin = sqrt(inverseOfTotalNumbersTested);
		double differenceFromExpectedPeriod2 = 
			fabs(averageAbsoluteDifferencePeriod2 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod4 = 
			fabs(averageAbsoluteDifferencePeriod4 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod8 = 
			fabs(averageAbsoluteDifferencePeriod8 - expectedAverageAbsoluteDifference);
		double differenceFromExpectedPeriod16 = 
			fabs(averageAbsoluteDifferencePeriod16 - expectedAverageAbsoluteDifference);

		result &= (differenceFromExpectedPeriod2 <= margin);
		result &= (differenceFromExpectedPeriod4 <= margin);
		result &= (differenceFromExpectedPeriod8 <= margin);
		result &= (differenceFromExpectedPeriod16 <= margin);

		//Log and return result:

		if (log) {
			printf("\nPRNG test (draw4spcg32s): %d drawn (normalized avg: %f)\n",
				                                               totalDrawn, avarege);
			printf("\tTiming: %f nanos per prn\n", nanosPerDraw);
			printf("\tZeroes: %d (expected: %f)\n", zeroes, expectedZeroes);
			printf("\tAvg. differences expected: %f (margin used: %f). Found:\n",
								expectedAverageAbsoluteDifference, margin);
			printf("\t\tPeriod 2: %f; Period 4: %f, Period 8: %f, Period 16: %f\n",
				averageAbsoluteDifferencePeriod2, averageAbsoluteDifferencePeriod4,
			   averageAbsoluteDifferencePeriod8, averageAbsoluteDifferencePeriod16);
		}

		if(!result) { LOG_ERROR("Something went wrong!"); }
		else if(log) { LOG_INFO("Passed"); }

		return result;
	}

	//used by the following test
	void makeParallelLoad(float loadLevel, std::chrono::microseconds testSleepTime);
	bool g_stopParallelLoad = false;
	inline bool sleepAndLog(std::chrono::microseconds sleepTestTime, int loadLevel, 
		              int thresholdLevel, std::chrono::microseconds* longestSnooze, 
                 std::chrono::microseconds accSnoozeTimes[], double baselineMicros,
  						    int* longestSnoozeLoadAndThresholdIndex, double margin);
						               
	#define SLEEP_TIME_MULT 4
	#define QUANTITY_LOADS 3
	#define QUANTITY_THRESHOLDS 3
	#define LOAD_THRESHOLD_COMBINATIONS (QUANTITY_LOADS*QUANTITY_THRESHOLDS)
	#define SLEEP_REPETITIONS 10
	#define RUNS_FOR_BASELINE 25
	#define BASELINE_MULTIPLIER_GRACE 10	

    //Sleeps from minimumSleepTimeMicros to maximumSleepTimeMicros, increasing by
	//multiplications of 4. Always includes maximumSleepTimeMicros.
	//For each sleep time, tests 3 different thresholds and 3 different parallel loads:
	//- No threshold, standard threshold and "full threshold" (busywaiting, as a baseline);
	//- No parrallel load, moderate load (~35% cpu time), high load (~90% cpu time);
	//For each test, calculates how much time was slept over what was expected.
	//Returns the proportion of tests which had excess <= margin*(sleeptime + baseline).
	//If log == true, logs to output (eg, to console) the results of each test.
	//Has default values for all parameters.
	static double testHybridBusySleeping(
			bool log = true, int minimumSleepTimeMicros = 17, 
		    int maximumSleepTimeMicros = 17408, double margin = 0.005) {

		LOG_DEBUG("Will test sleeping functionality", 1);

		float loads[QUANTITY_LOADS];
		loads[0] = 0;
		loads[1] = 0.35f;
		loads[2] = 0.9f;

		double maxMinMultiplier = maximumSleepTimeMicros / minimumSleepTimeMicros;
		int sleepSteps = (int)ceil(log2(maxMinMultiplier)/2);

		std::thread* parallelLoad_ptr;
		
		std::chrono::microseconds minSleepTime = 
				std::chrono::microseconds(minimumSleepTimeMicros);
		std::chrono::microseconds maxSleepTime = 
				std::chrono::microseconds(maximumSleepTimeMicros);
		std::chrono::microseconds timeToSleep;

		std::chrono::microseconds accSnoozeTimes[LOAD_THRESHOLD_COMBINATIONS];
		for (int i = 0; i < LOAD_THRESHOLD_COMBINATIONS; i++) {
			accSnoozeTimes[i] = std::chrono::microseconds(0);
		}
		std::chrono::microseconds longestSnooze = std::chrono::microseconds(0);
		int longestSnoozeLoadAndThresholdIndex = -1;

		//int sleptMoreThanExpected = 0; //TODO: Review and delete

		std::chrono::microseconds start;
		std::chrono::microseconds end;
		std::chrono::microseconds snooze;
		double baselineMicros = 0;
		for (int i = 0; i < RUNS_FOR_BASELINE; i++) {
			start = nowMicros();
			hybridBusySleepForMicros(minSleepTime, 2*minSleepTime);
			end = nowMicros();
			snooze = (end - start) - minSleepTime;
			baselineMicros += (double)snooze.count()/RUNS_FOR_BASELINE;
		}

		int passed = 0;

		for (int i = 0; i < sleepSteps; i++) {
			for (int j = 0; j < QUANTITY_LOADS; j++) {

				timeToSleep = minSleepTime*(int)round(pow(SLEEP_TIME_MULT, i));
				if (timeToSleep > maxSleepTime) {
					timeToSleep = maxSleepTime;
				}

				g_stopParallelLoad = false;
				parallelLoad_ptr = new std::thread(makeParallelLoad, loads[j], timeToSleep);
				
				for(int k = 0; k < SLEEP_REPETITIONS; k++){

					passed += (int)sleepAndLog(timeToSleep, j, 0, &longestSnooze, 
						                          accSnoozeTimes, baselineMicros, 
						             &longestSnoozeLoadAndThresholdIndex, margin);

					passed += (int)sleepAndLog(timeToSleep, j, 1, &longestSnooze, 
						                          accSnoozeTimes, baselineMicros, 
						             &longestSnoozeLoadAndThresholdIndex, margin);

					passed += (int)sleepAndLog(timeToSleep, j, 2, &longestSnooze, 
						                          accSnoozeTimes, baselineMicros, 
						             &longestSnoozeLoadAndThresholdIndex, margin);
				}
				
				g_stopParallelLoad = true;
				parallelLoad_ptr->join();
			}
		}

		int totalCycles = QUANTITY_LOADS*QUANTITY_THRESHOLDS*SLEEP_REPETITIONS*sleepSteps;
		double proportionPassed = (double)passed/totalCycles;
		int longestSnoozeThresholdLvl = longestSnoozeLoadAndThresholdIndex % QUANTITY_THRESHOLDS;
		int longestSnoozeParallelLoadLvl = longestSnoozeLoadAndThresholdIndex / QUANTITY_THRESHOLDS;
		
		double totalAvgSnooze = 0;
		double avgSnoozeTimes[LOAD_THRESHOLD_COMBINATIONS];
		for (int i = 0; i < LOAD_THRESHOLD_COMBINATIONS; i++) {
			avgSnoozeTimes[i] = (double)accSnoozeTimes[i].count()/(SLEEP_REPETITIONS*sleepSteps);
			totalAvgSnooze += avgSnoozeTimes[i]/LOAD_THRESHOLD_COMBINATIONS;
		}		

		if (log) {
			printf("\nHybrid Sleep test completed. Ran %d total cycles (Passed %f%%), target sleep: from %d to %d microseconds\n",
					totalCycles, proportionPassed*100, minimumSleepTimeMicros, maximumSleepTimeMicros);
			printf("\tAvg snooze: %f microseconds (baseline: %f microseconds). Longest: %ld microseconds (thresh lvl %d, load lvl %d)\n",
					totalAvgSnooze, baselineMicros, longestSnooze.count(),
				    longestSnoozeThresholdLvl, longestSnoozeParallelLoadLvl);
			for (int i = 0; i < QUANTITY_LOADS; i++) {
				for (int j = 0; j < QUANTITY_THRESHOLDS; j++) {
					printf("\t\tFor load lvl %d, threshold level %d: avg snooze: %f microseconds\n",
							i, j, avgSnoozeTimes[QUANTITY_THRESHOLDS*i + j]);
				}
			}

		}
		
		return proportionPassed;
	}

	void makeParallelLoad(float loadLevel, std::chrono::microseconds testSleepTime) {

		float cyclesPerTestSleepTime = 1/loadLevel; //more load, more cycles

		int sleepMicros = (int)round((1 - loadLevel)*testSleepTime.count()/cyclesPerTestSleepTime);
		int spinMicros = (int)round(loadLevel*testSleepTime.count()/cyclesPerTestSleepTime);		

		std::chrono::microseconds sleepTime = std::chrono::microseconds(sleepMicros);
		std::chrono::microseconds spinTime = std::chrono::microseconds(spinMicros);

		while (!g_stopParallelLoad) {
			hybridBusySleepForMicros(sleepTime);
			hybridBusySleepForMicros(spinTime, 2*spinTime); //busy spins from beggining
		}
	}

	inline bool sleepAndLog(std::chrono::microseconds sleepTestTime, int loadLevel, 
		              int thresholdLevel, std::chrono::microseconds* longestSnooze, 
				 std::chrono::microseconds accSnoozeTimes[], double baselineMicros, 
  						    int* longestSnoozeLoadAndThresholdIndex, double margin){

		std::chrono::microseconds start;
		std::chrono::microseconds end;

		if(thresholdLevel == 0){
			start = nowMicros();
			hybridBusySleepForMicros(sleepTestTime, std::chrono::microseconds(0));
			end = nowMicros();
		}
		else if(thresholdLevel == 1){
			start = nowMicros();
			hybridBusySleepForMicros(sleepTestTime);
			end = nowMicros();
		}
		else if(thresholdLevel == 2){
			start = nowMicros();
			hybridBusySleepForMicros(sleepTestTime, 2*sleepTestTime); //busy spins from beggining
			end = nowMicros();
		}
		else{
			return false;
		}
		
		auto snooze = (end - start) - sleepTestTime;
		int loadAndThresholdIndex = QUANTITY_THRESHOLDS*loadLevel + thresholdLevel;
		accSnoozeTimes[loadAndThresholdIndex] += (snooze);
		if (snooze > *longestSnooze) {
			*longestSnooze = snooze;
			*longestSnoozeLoadAndThresholdIndex = loadAndThresholdIndex;
		}

		bool result = true;
		if ((double)snooze.count() > (BASELINE_MULTIPLIER_GRACE * baselineMicros)) {
			double effectiveSnoozeOverBaseline = (double)snooze.count() - baselineMicros;
			result = effectiveSnoozeOverBaseline <= (margin*(double)sleepTestTime.count());
		}

		return result;
	}

	#define MINIMUM_PROPORTION_SLEEP_PASSES (0.95)
	static bool testSnooze(bool printLog = true) {
	
		LOG_DEBUG("Will test sleeping and waking a few times...\n", 1); 
		GETCHAR_PAUSE; 

		double result = AZ::testHybridBusySleeping(printLog); 

		bool passed = result > MINIMUM_PROPORTION_SLEEP_PASSES;
		if(!passed) { LOG_ERROR("Snoozed more than the maximum margin set"); }
		else if(printLog) { LOG_INFO("Passed"); }

		return passed;
	}
}
