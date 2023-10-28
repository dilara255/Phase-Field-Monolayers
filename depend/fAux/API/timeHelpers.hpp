#pragma once

#include "miscStdHeaders.h"

#define AZ_MICROS_TO_BUSY_WAIT 25

namespace AZ{

	const static std::chrono::nanoseconds zeroNano = std::chrono::nanoseconds(0);
	const static std::chrono::nanoseconds zeroMicro = std::chrono::microseconds(0);
	const static std::chrono::nanoseconds zeroMilli = std::chrono::milliseconds(0);
	const static std::chrono::nanoseconds zeroSeconds = std::chrono::seconds(0);

	const static std::chrono::microseconds standardThreshold(AZ_MICROS_TO_BUSY_WAIT);

	//takes about 0.22 microseconds per call on test system
	inline std::chrono::microseconds nowMicros() {
		auto now = std::chrono::steady_clock::now().time_since_epoch();
		return std::chrono::duration_cast<std::chrono::microseconds>(now);
	}

	inline std::chrono::nanoseconds nowNanos() {
		auto now = std::chrono::steady_clock::now().time_since_epoch();
		return std::chrono::duration_cast<std::chrono::nanoseconds>(now);
	}

	//This loops while yielding via Sleep(0) up to targetWakeTime - threshold,
	//then busy-waits until targetWakeTime. Does nothing if now >= targetWakeTime.
	//If no threshold is especified, uses standard value defined by Aux0.
	//WARNING: NO IDEA how portable this is (is it HACKY?)
	//TODO: test response for different load scenarios
	//TODO: run tests on other platforms
	void hybridBusySleep(std::chrono::steady_clock::time_point targetWakeTime,
		                 std::chrono::microseconds threshold = standardThreshold);
	//This loops while yielding via Sleep(0) for microsToSleep - threshold,
	//then busy-waits to complete microsToSleep. 
	//If no threshold is especified, uses standard value defined by Aux0.
	void hybridBusySleepForMicros(std::chrono::microseconds microsToSleep,
		                          std::chrono::microseconds threshold = standardThreshold);
}