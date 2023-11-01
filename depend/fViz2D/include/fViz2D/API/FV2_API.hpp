#pragma once

#include "fAux/API/core.hpp"

#include "fViz2D/utils/imageUtils.hpp"
#include "fViz2D/utils/guiUtils.hpp"
#include "fViz2D/renderer/rendererData.hpp"

#include "fViz2D/resourcePaths.hpp"

namespace F_V2 {
	F_V2_API std::thread spawnRendererOnNewThread(bool* externalBool_ptr, IMG::rgbaImage_t* dynamicData_ptr, 
									          COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr, 
		                                      F_V2::rendererRetCode_st* returnCode_ptr, 
		                                      const char* bannerPathFromBinary  = F_V2::testBannerPathFromBinary);

	//TODO: ADD TEST
	F_V2_API F_V2::rendererRetCode_st spawnRendererOnThisThread(bool* externalBool_ptr, 
		                               IMG::rgbaImage_t* dynamicData_ptr, COLOR::rgbaF_t* clearColor_ptr, 
		                               COLOR::rgbaF_t* noiseTint_ptr, 
		                               const char* bannerPathFromBinary  = F_V2::testBannerPathFromBinary);

	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//TODO: CRITICAL: WARNING: THIS-IS-SO-UGLY :' ):
	//load4channelTextureFromRgbaImage vs loadR32FtextureFromFloats need of IMG::rgbaImage_t* vs IMG::floats2Dfield_t*
	//is the only thing keeping me from having all of this be a single thing, receiving a more general "IMG::2Ddata_t*"
	//MAKE THAT WORK BEFOR BUILDING NEW FEATURES, PLS
	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	F_V2_API std::thread spawnRendererOnNewThreadF(bool* externalBool_ptr, IMG::floats2Dfield_t* dynamicData_ptr, 
									          COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr, 
		                                      F_V2::rendererRetCode_st* returnCode_ptr, 
		                                      const char* bannerPathFromBinary  = F_V2::testBannerPathFromBinary);

	//TODO: ADD TEST
	F_V2_API F_V2::rendererRetCode_st spawnRendererOnThisThreadF(bool* externalBool_ptr, 
		                               IMG::floats2Dfield_t* dynamicData_ptr, COLOR::rgbaF_t* clearColor_ptr, 
		                               COLOR::rgbaF_t* noiseTint_ptr, 
		                               const char* bannerPathFromBinary  = F_V2::testBannerPathFromBinary);
}