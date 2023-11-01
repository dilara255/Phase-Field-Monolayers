#pragma once

namespace F_V2 {
	enum struct rendererRetCode_st {OK, STILL_RUNNING, CONTEXT_ACQ_FAILED, IMAGE_LOAD_FAILED, DYNAMIC_IMAGE_INITIAL_LOAD_FAILED};
	enum struct texRetCode_st {OK, SIZES_DONT_MATCH_FOR_COPY};
}