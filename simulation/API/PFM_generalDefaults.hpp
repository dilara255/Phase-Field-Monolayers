#pragma once

//Default values which *do not* depend on PFM specific data definitions
//TODO: brief summary

#include "PFM_API_enums.hpp"

#include "fAux/API/miscStdHeaders.h"

namespace PFM {

	static const char* argumentNames[TOTAL_MAIN_ARGS] = { "ProgramCall", "SimFuncToRun(uint)", "Lamba(double)", 
		                                      "Gamma(double)", "dt(double)", "Cells(uint)", "Width(uint)", "Height(uint)", 
											  "InitialCondition(uint)", "Bias(double)", "Seed(uint64)", 
											  "Method(uint)", "StartPaused(bool)", 
											  "ChangePerElementPerStepToStop(double)", "MaximumSteps(uint64)",
	                                          "StepsPerCheck(uint)", "ChangePerCheck(double)",
		                                      "CallerKey(uint)", "AdaptativeDt(bool)", "MaxAvgChangePerStep(double)",
											  "MaxSpeedupMult(double)", "MinSlowdownMult(double)", 
											  "UseMaxSafeDt(bool)"};
	static const char* defaultArgument = "default";

	static const double lowStableEq = 0;
	static const double highStableEq = 1;

	static const uint64_t completelyArbitraryStepToUnlockFullDt = 1000;

	//For adaptative dt:
	
	static const double distanceStableEquilibria = highStableEq - lowStableEq;
	static const double defaultMaxChangePerStep = PFM::distanceStableEquilibria / 8.0;
	static const double defaultMaxSpeedupMult = 1.00001; //note that larger values may make the average *slower*
	static const double defaultMinSlowdownMult = 0.97; //because it might cause more drastic bump downs (to avoid crashes)

	//"Margin" from the equilibrium points away from the interface. 
	//"Climbing" the potential too far outside from this may cause numerical instability.
	static const double outsideLimits = distanceStableEquilibria/10;

	//Maps [-this, 1 + this] to PRGMs [0,255], to see "around" the wells
	static const double defaultPGMmargin = outsideLimits;

	static const uint32_t defaulStepsPerCheck = 5000;
	static const double defaultAbsChangePerCheck = 0.0025;
	
	static const double defaultAbsChangePerStepToStop = 0.0000000015;
	static const uint64_t defaulMaxSteps = 6000000;
	static const int stepsAtChangeThresholdToStop = 100; 

	static const uint32_t defaultCallerKey = 0;
}