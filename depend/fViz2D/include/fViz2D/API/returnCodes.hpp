#pragma once

namespace F_V2 {
	enum struct rendererRetCode_st { OK, STILL_RUNNING, CONTEXT_ACQ_FAILED, IMAGE_LOAD_FAILED, 
		                            DYNAMIC_IMAGE_INITIAL_LOAD_FAILED, BAD_DYNAMIC_DATA_FORMAT };
	enum struct texRetCode_st { OK, SIZES_DONT_MATCH_FOR_COPY, BAD_VALUE_FIELD };
	enum struct imageFileRetCode_st { OK, SAVE_FAILED, BAD_DATA_TO_SAVE };
}