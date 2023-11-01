#include "fAux/API/miscStdHeaders.h"
#include "fAux/API/logAPI.hpp"
#include "fAux/API/prng.hpp"
#include "fAux/API/timeHelpers.hpp"

#include "fViz2D/API/FV2_testsAPI.hpp"
#include "fViz2D/API/FV2_API.hpp"

#include "PFM_API.hpp"

#include "guiTests.hpp"

bool PFM_GUI_TESTS::runFviz2DinternalTests() {
	
	bool result = F_V2::linkAndLogTest();
	result &= F_V2::rendererTestFromImage();
	return result & F_V2::rendererTestFromDoubles2Dfield();
}

bool PFM_GUI_TESTS::fViz2DintegrationTests() {
	LOG_DEBUG("PFM_GUI <-> fViz integration tests");
	LOG_DEBUG("Same thread test (red static noise expected");

	IMG::floats2Dfield_t field = IMG::createFloats2Dfield(512, 512);
	uint64_t seed = DEFAULT_PRNG_SEED0;
	for (size_t i = 0; i < field.size.getTotalElements(); i++) {
		field.data[i] = (float)(AZ::draw1spcg32(&seed)/(double)UINT32_MAX);
	}

	COLOR::rgbaF_t clear = COLOR::CLEAR;
	COLOR::rgbaF_t tint = COLOR::FULL_WHITE;

	bool works = false;
	F_V2::rendererRetCode_st retCode = F_V2::spawnRendererOnThisThreadF(&works, &field, &clear, &tint);

	bool result1 = (retCode == F_V2::rendererRetCode_st::OK && works);
	if (result1) { LOG_INFO("OK"); }
	else { LOG_ERROR("ERROR"); }

	GETCHAR_PAUSE;

	LOG_DEBUG("Multi-thread test (red dynamic noise expected");

	retCode = F_V2::rendererRetCode_st::STILL_RUNNING;
	works = false;
	std::thread renderTrhread = F_V2::spawnRendererOnNewThreadF(&works, &field, &clear, &tint, &retCode);

	while (retCode == F_V2::rendererRetCode_st::STILL_RUNNING) {
		for (size_t i = 0; i < field.size.getTotalElements(); i++) {
			field.data[i] += 0.01 * (float)(AZ::draw1spcg32(&seed)/(double)UINT32_MAX);
			field.data[i] -= 1 * (field.data[i] > 1.0);
		}

		AZ::hybridBusySleepForMicros(std::chrono::microseconds(1000));
	}

	renderTrhread.join();

	bool result2 = (retCode == F_V2::rendererRetCode_st::OK && works);
	if (result2) { LOG_INFO("OK"); }
	else { LOG_ERROR("ERROR"); }

	GETCHAR_PAUSE;

	return result1 & result2;
}

bool PFM_GUI_TESTS::guiLinkingAndDependencyTests() {
	PFM::linkingTest();	
	
	bool result1 = PFM_GUI_TESTS::runFviz2DinternalTests();
	if (result1) { LOG_INFO("fViz2D internal tests passed"); }
	else { LOG_ERROR("fViz2D internal tests NOT passed"); }
	
	bool result2 = PFM_GUI_TESTS::fViz2DintegrationTests();
	if (result2) { LOG_INFO("fViz2D integration tests passed"); }
	else { LOG_ERROR("fViz2D integration tests NOT passed"); }

	if (result1 & result2) { LOG_INFO("All tests passed!"); }
	else { LOG_ERROR("Some of the tests failed"); }

	GETCHAR_FORCE_PAUSE;

	return (result1 & result2);
}