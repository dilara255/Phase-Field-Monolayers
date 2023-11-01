#pragma once

//TODO: review: some stuff here probably could be made more general to support more types (template?)

#include "GLFW/glfw3.h"
#include "stbImage/stb_image.h"

#include "fViz2D/API/returnCodes.hpp"

namespace IMG {

	typedef struct imgSizeInfo_st {
		size_t width = 0, height = 0;
		size_t channels;
		size_t bytesPerChannel;
		bool initialized = false;

		size_t getStride() { return channels * width; }
		size_t getIndex(uint32_t row, uint32_t collumn, uint8_t channel) {
			return (row * getStride()) + (collumn * channels) + channel;
		}
		size_t getMaxIndex() { return ( (channels * width * height) - 1 ); }
		size_t getTotalElements() { return channels * width * height; }
	} imgSizeInfo_t;

	typedef struct rgbaImage_st {
		unsigned char* data;
		imgSizeInfo_t size = {0, 0, 4, 1};
		bool initialized = false;

		~rgbaImage_st() { stbi_image_free(data); }
	} rgbaImage_t;

	// Loads rgbaImage_t from a file. The data buffer becomes a responsability of the caller.
	rgbaImage_t load4channel8bpcImageFromFile(const char* filename);

	typedef struct doubles2Dfield_st {
		double* data;
		imgSizeInfo_t size = {0, 0, 1, 8};

		~doubles2Dfield_st() { free(data); data = NULL; }
	} doubles2Dfield_t;

	typedef struct floats2Dfield_st {
		float* data;
		imgSizeInfo_t size = {0, 0, 1, 4};

		~floats2Dfield_st() { free(data); data = NULL; }
	} floats2Dfield_t;

	//In case alocation fails, the field's "initialized" member will be false and width/height will be zero
	doubles2Dfield_t createDoubles2Dfield(size_t width, size_t height);
	floats2Dfield_t createFloats2Dfield(size_t width, size_t height);

	//Expects origin and destination to have the same amount of elements
	F_V2::texRetCode_st copy2Dfield(const floats2Dfield_t* floatOrigin_ptr, doubles2Dfield_t* doubleDest_ptr);
	F_V2::texRetCode_st copy2Dfield(const doubles2Dfield_t* doubleOrigin_ptr, floats2Dfield_t* floatDest_ptr);
}

namespace COLOR {
	
	typedef enum {R, G, B, A} channel_t;

	typedef struct rgbF_st {
		float r, g, b;

		float getChannel(channel_t channel) {
			switch (channel) {
				case channel_t::R: return r;
				case channel_t::G: return g;
				case channel_t::B: return b;
				case channel_t::A: return (r *g * b); //TODO: this is ugly
			}
		}
	} rgbF_t;

	typedef struct rgbaF_st {
		float r, g, b, a;

		float getChannel(channel_t channel) {
			switch (channel) {
				case channel_t::R: return r;
				case channel_t::G: return g;
				case channel_t::B: return b;
				case channel_t::A: return a;
			}
		}
	} rgbaF_t;

	typedef struct rgbaC_st {
		char r, g, b, a;

		char getChannel(channel_t channel) {
			switch (channel) {
				case channel_t::R: return r;
				case channel_t::G: return g;
				case channel_t::B: return b;
				case channel_t::A: return a;
			}
		}
	} rgbaC_t;



	//TODO: a reasonable way to have defaults for both formats
	const rgbaF_t CLEAR = {0.263f, 0.141f, 0.384f, 1.0f};
	const rgbaF_t FULL_WHITE = {1.0f, 1.0f, 1.0f, 1.0f};
}
