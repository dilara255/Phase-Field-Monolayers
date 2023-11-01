#pragma once

#include "utils/imageUtils.hpp" //TODO: I don't want this to include any non-api header
#include "returnCodes.hpp"

namespace F_V2 {
	F_V2::rendererRetCode_st rendererMain(bool* externalBool_ptr, IMG::rgbaImage_t* dynamicData_ptr, 
		             COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr,
		             const char* bannerPathFromBinary);

	void rendererMainForSeparateThread(bool* externalBool_ptr, IMG::rgbaImage_t* dynamicData_ptr, 
									   COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr, 
		                               F_V2::rendererRetCode_st* returnCode_ptr, 
		                               const char* bannerPathFromBinary);

	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//TODO: CRITICAL: WARNING: THIS-IS-SO-UGLY :' ):
	//load4channelTextureFromRgbaImage vs loadR32FtextureFromFloats need of IMG::rgbaImage_t* vs IMG::floats2Dfield_t*
	//is the only thing keeping me from having all of this be a single thing, receiving a more general "IMG::2Ddata_t*"
	//MAKE THAT WORK BEFOR BUILDING NEW FEATURES, PLS
	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	F_V2::rendererRetCode_st rendererMainF(bool* externalBool_ptr, IMG::floats2Dfield_t* dynamicData_ptr, 
		             COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr,
		             const char* bannerPathFromBinary);

	void rendererMainForSeparateThreadF(bool* externalBool_ptr, IMG::floats2Dfield_t* dynamicData_ptr, 
									   COLOR::rgbaF_t* clearColor_ptr, COLOR::rgbaF_t* noiseTint_ptr, 
		                               F_V2::rendererRetCode_st* returnCode_ptr, 
		                               const char* bannerPathFromBinary);
}