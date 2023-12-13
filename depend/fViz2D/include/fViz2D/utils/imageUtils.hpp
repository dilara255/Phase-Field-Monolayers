#pragma once

//TODO: review: some stuff here probably could be made more general to support more types (template?)
//TODO-CRITICAL: write tests for this stuff
//TODO: maybe generic data type could be generic only on the pointer to data, and hold size directly?
//TODO: some stuff here should eventually be templated. Also take a look at the alternative architecture (loose file)

#include "GLFW/glfw3.h"
#include "stbImage/stb_image.h"
#include "stbImage/stb_image_write.h"

#include "fAux/API/miscStdHeaders.h"

#include "fViz2D/API/returnCodes.hpp"

#include <assert.h>
#include <memory>
#include <vector>

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
		unsigned char r, g, b, a;

		char getChannel(channel_t channel) {
			switch (channel) {
				case channel_t::R: return r;
				case channel_t::G: return g;
				case channel_t::B: return b;
				case channel_t::A: return a;
			}
		}
	} rgbaC_t;

	typedef unsigned char grey8bpcC_t;

	//TODO: a reasonable way to have defaults for both formats
	const rgbaF_t CLEAR = {0.263f, 0.141f, 0.384f, 1.0f};
	const rgbaF_t FULL_WHITE = {1.0f, 1.0f, 1.0f, 1.0f};
	const rgbaF_t FULL_BLACK = {0.0f, 0.0f, 0.0f, 1.0f};
	const rgbaF_t DEBUG_PINK = {1.0f, (1.0f/3), (13.0f/15), 1.0f};

	const rgbaC_t CLEAR_8B = {67, 36, 98, 255};
	const rgbaC_t FULL_WHITE_8B = {255, 255, 255, 255};
	const rgbaC_t FULL_BLACK_8B = {0, 0, 0, 255};
	const rgbaC_t FULL_BLUE_8B = {0, 0, 255, 255};
	const rgbaC_t DARK_BLUE_8B = {0, 0, 127, 255};
	const rgbaC_t FULL_YELLOW_8B = {255, 255, 0, 255};
	const rgbaC_t FULL_GREEN_8B = {0, 255, 0, 255};
	const rgbaC_t FULL_RED_8B = {255, 0, 0, 255};
	const rgbaC_t DARK_RED_8B = {127, 0, 0, 255};
	const rgbaC_t DEBUG_PINK_8B = {255, 85, 221, 255};

	const grey8bpcC_t DEBUG_GRAY = 42;

	typedef struct valueRgba8bpcColorPair_st {
		double value;
		COLOR::rgbaC_t color;
	} valueRgba8bpcColorPair_t;

	typedef std::vector<valueRgba8bpcColorPair_t> schemeVectorRgba8bpc_t;
	
	static const schemeVectorRgba8bpc_t defaultPassThroughScheme = {};
	static const schemeVectorRgba8bpc_t defaultBWscheme = { valueRgba8bpcColorPair_t{ 0.0, FULL_BLACK_8B }, 
		                                            valueRgba8bpcColorPair_t{ 1.0, FULL_WHITE_8B } };
	static const schemeVectorRgba8bpc_t defaultBlueRedScheme = { valueRgba8bpcColorPair_t{-0.5, FULL_BLUE_8B }, 
		                                                 valueRgba8bpcColorPair_t{ 1.5, FULL_RED_8B } };
	static const schemeVectorRgba8bpc_t defaultBlueYellowRedScheme = { valueRgba8bpcColorPair_t{-0.5, FULL_BLACK_8B }, 
															   valueRgba8bpcColorPair_t{ 0.0, FULL_BLUE_8B },  
		                                                       valueRgba8bpcColorPair_t{ 0.5, FULL_YELLOW_8B }, 
															   valueRgba8bpcColorPair_t{ 1.0, FULL_RED_8B }, 
		                                                       valueRgba8bpcColorPair_t{ 1.5, FULL_WHITE_8B } };

	//TODO: Add schemes for other formats

	//TODO: actual support for shader schemes not yet implemented
	typedef struct schemeRgba8bpcForShader_st {
		int totalElements = 0;
		std::vector<float> values;
		std::vector<COLOR::rgbaC_t> colors;
	} schemeRgba8bpcForShader_t;

	//Insertions guarantee that values are in increasing order
	//TODO: this has grown and should be a class
	//TODO-CRITICAL: also please move the definitions to the .cpp file
	typedef struct colorInterpolation_st {
		private:
			schemeVectorRgba8bpc_t m_correspondences;
			double m_span = 0;
			double m_bias = 0;
		
	public:
			const schemeVectorRgba8bpc_t* const correspondences_ptr = &m_correspondences;
			const double* const bias_ptr = &m_bias;
			const double* const span_ptr = &m_span;

			bool addCorrespondence(double newValue, COLOR::rgbaC_t newColor) {  
				if(m_correspondences.size() == 0 || m_correspondences.back().value < newValue) { 
					m_correspondences.push_back({newValue, newColor});
					m_bias = m_correspondences.at(0).value;
					m_span = newValue - m_bias;
					
					return true; 
				}
				else { return false; }
			}

			void loadScheme(const schemeVectorRgba8bpc_t* schemeToLoad_ptr) { 
				if(!m_correspondences.empty()) { m_correspondences.clear(); }
				m_correspondences = *schemeToLoad_ptr;
				if(!m_correspondences.empty()) {
					m_bias = m_correspondences.at(0).value;
					m_span = m_correspondences.back().value - m_bias;
				}
				else { m_bias = 0.0, m_span = 0.0; }
			}

			void changeBiasTo(double newBias) {
				double change = newBias - m_bias;
				for(size_t i = 0; i < m_correspondences.size(); i++) { m_correspondences.at(i).value += change; }
				m_bias = newBias;
			}

			//Must be > 0 and only works if the scheme already has two or more elements
			void changeSpanTo(double newSpan) {				
				if(newSpan <= 0 || m_correspondences.size() < 2) { return; }

				assert(m_span > 0);

				double change = newSpan/m_span;
				for(size_t i = 0; i < m_correspondences.size(); i++) { 
					double oldValue = m_correspondences.at(i).value;
					m_correspondences.at(i).value = m_bias + change*(oldValue - m_bias); 
				}
				m_span = newSpan;
			}

			//Sets the values interval to [0.0, 1.0]
			void normalizeSpan() { changeBiasTo(0.0); changeSpanTo(1.0); }

			//Sets up the schemeRgba8bpcForShader_t pointed at by shaderScheme_ptr as a normalized copy of this
			//Deletes any old data on shaderScheme_ptr
			void setSchemeForShader(schemeRgba8bpcForShader_t* shaderScheme_ptr) {
				
				shaderScheme_ptr->values.clear();
				shaderScheme_ptr->colors.clear();

				shaderScheme_ptr->totalElements = m_correspondences.size();
				shaderScheme_ptr->values.reserve(shaderScheme_ptr->totalElements);
				shaderScheme_ptr->colors.reserve(shaderScheme_ptr->totalElements);

				colorInterpolation_st tempScheme = *this;
				tempScheme.normalizeSpan();

				for (int i = 0; i < shaderScheme_ptr->totalElements; i++) {
					shaderScheme_ptr->values.at(i) = (float)tempScheme.m_correspondences.at(i).value;
					shaderScheme_ptr->colors.at(i) = tempScheme.m_correspondences.at(i).color;
				}
			}

	} colorInterpolation_t;
		
	rgbaC_t interpolateTwoColors(double t, const rgbaC_t* colorBefore_ptr, const rgbaC_t* colorAfter_ptr);

	//Returns DEBUG_PINK_8B if an empty scheme is sent
	rgbaC_t interpolatedColorFromValue(double value, const  colorInterpolation_t* scheme_ptr);
}

namespace IMG {

	typedef struct imgSizeInfo_st {
		size_t width = 0, height = 0;
		size_t channels;
		size_t bytesPerChannel;
		bool initialized = false;

		size_t getStride() const { return channels * width; }
		size_t getStrideInBytes() const { return channels * width * bytesPerChannel; }
		//Will return a bad index if bad parameters are passed on release
		size_t getIndex(uint32_t row, uint32_t collumn, uint8_t channel) const {
			assert( (row < height) && (collumn < width) && (channel < channels) );
			return (row * getStride()) + (collumn * channels) + channel;
		}
		//Will return a bad index if bad parameters are passed on release
		size_t getLinearIndexOfChannel(size_t element, size_t channel) const {
			assert( (element < height * width) && (channel < channels) );
			return ( (channels * element) + channel);
		}
		size_t getMaxIndex() const { return ( (channels * width * height) - 1 ); }
		size_t getTotalElements() const { return channels * width * height; }
		size_t getTotalArea() const { return width * height; }
	} imgSizeInfo_t;

	typedef struct grey8bpcImage_st {
		std::unique_ptr<unsigned char> data = NULL;
		imgSizeInfo_t size = {0, 0, 1, 1};
	} grey8bpcImage_t;

	//In case alocation fails, the field's "initialized" member will be false and width/height will be zero
	grey8bpcImage_t createEmpty1channel8bpcImage(size_t width, size_t height);

	typedef struct rgbaImage_st {
		std::unique_ptr<unsigned char> data = NULL;
		imgSizeInfo_t size = {0, 0, 4, 1};
	} rgbaImage_t;

	// Loads rgbaImage_t from a file. The data buffer becomes a responsability of the caller.
	rgbaImage_t load4channel8bpcImageFromFile(const char* filename);
	//In case alocation fails, the field's "initialized" member will be false and width/height will be zero
	rgbaImage_t createEmpty4channel8bpcImage(size_t width, size_t height);

	typedef struct floats2Dfield_st {
		std::unique_ptr<float> data = NULL;
		imgSizeInfo_t size = {0, 0, 1, 4};
	} floats2Dfield_t;

	typedef struct doubles2Dfield_st {
		std::unique_ptr<double> data = NULL;
		imgSizeInfo_t size = {0, 0, 1, 8};
	} doubles2Dfield_t;

	//TODO: using union to represent general data type. Possibly should turn these into classes
	enum class kinds2Ddata { UNINITIALIZED_UNION, GREY_8BPC, RGBA_IMAGE, FLOATS_FIELD, DOUBLES_FIELD };
	union generic2DfieldPtrs_u { grey8bpcImage_t* grey8bpc_ptr = nullptr; rgbaImage_t* rgbaField_ptr; 
	                             floats2Dfield_t* floatsField_ptr; doubles2Dfield_t* doublesField_ptr; };

	typedef struct generic2DfieldPtr_st {

		void storeGrey8bpc(grey8bpcImage_t* grey8bpc_ptr) { 
			field_ptr.grey8bpc_ptr = grey8bpc_ptr; 
			kindOfField = kinds2Ddata::GREY_8BPC; 
		}

		void storeRGBAfield(rgbaImage_t* rgbaField_ptr) { 
			field_ptr.rgbaField_ptr = rgbaField_ptr; 
			kindOfField = kinds2Ddata::RGBA_IMAGE; 
		}

		void storeFloatsField(floats2Dfield_t* floatsField_ptr) { 
			field_ptr.floatsField_ptr = floatsField_ptr; 
			kindOfField = kinds2Ddata::FLOATS_FIELD; 
		}

		void storeDoublesField(doubles2Dfield_t* doublesField_ptr) { 
			field_ptr.doublesField_ptr = doublesField_ptr; 
			kindOfField = kinds2Ddata::DOUBLES_FIELD; 
		}

		generic2DfieldPtrs_u getFieldPtr() { return field_ptr; }
		const generic2DfieldPtrs_u getConstFieldPtr() const { return field_ptr; }
		kinds2Ddata getKindOfField() const { return kindOfField; }
		
		//Returns nullptr in case the field is not initialized
		const imgSizeInfo_t* getSizeInfo_ptr() const { 

			switch(kindOfField) {
				case kinds2Ddata::UNINITIALIZED_UNION: return nullptr;
				case kinds2Ddata::GREY_8BPC: return &(field_ptr.grey8bpc_ptr->size);
				case kinds2Ddata::RGBA_IMAGE: return &(field_ptr.rgbaField_ptr->size);
				case kinds2Ddata::FLOATS_FIELD: return &(field_ptr.floatsField_ptr->size);
				case kinds2Ddata::DOUBLES_FIELD: return &(field_ptr.doublesField_ptr->size);
				default:
					assert(false);
					return nullptr;
			} 
		}

		//Returns nullptr in case the field is not initialized
		const void* getVoidData_ptr() const { 

			switch(kindOfField) {
				case kinds2Ddata::UNINITIALIZED_UNION: return nullptr;
				case kinds2Ddata::GREY_8BPC: return (void*)(field_ptr.grey8bpc_ptr->data.get());
				case kinds2Ddata::RGBA_IMAGE: return (void*)(field_ptr.rgbaField_ptr->data.get());
				case kinds2Ddata::FLOATS_FIELD: return (void*)(field_ptr.floatsField_ptr->data.get());
				case kinds2Ddata::DOUBLES_FIELD: return (void*)(field_ptr.doublesField_ptr->data.get());
				default:
					assert(false);
					return nullptr;
			} 
		}

		private:
			kinds2Ddata kindOfField = kinds2Ddata::UNINITIALIZED_UNION;
			generic2DfieldPtrs_u field_ptr;

	} generic2DfieldPtr_t;
	

	//In case alocation fails, the field's "initialized" member will be false and width/height will be zero
	doubles2Dfield_t createDoubles2Dfield(size_t width, size_t height);
	floats2Dfield_t createFloats2Dfield(size_t width, size_t height);

	//Expects origin and destination to have the same amount of elements
	F_V2::texRetCode_st copy2Dfield(const floats2Dfield_t* floatOrigin_ptr, doubles2Dfield_t* doubleDest_ptr);
	F_V2::texRetCode_st copy2Dfield(const doubles2Dfield_t* doubleOrigin_ptr, floats2Dfield_t* floatDest_ptr);

	//Copies values to a greyscale image, preserving both. Maps values linearly to the interval [min, min+span]
	//Will return an empty grey8bpcImage_t if span <= 0, the origin is not initialized or something fails
	grey8bpcImage_t copy2DfieldToNewGreyscaleImage(const floats2Dfield_t* floatOrigin_ptr, double min, double span);
	grey8bpcImage_t copy2DfieldToNewGreyscaleImage(const doubles2Dfield_t* doubleOrigin_ptr, double min, double span);

	//These expects values and image fields to have the same amount of "pixels" (same width * height)
	//They also expects the scheme to have at least two points. Return OK on success
	//The generic versions considers it an error if valuesField_ptr is already an image type (rgba, greyscale, etc)
	F_V2::texRetCode_st translateValuesToInterpolatedColors(const generic2DfieldPtr_t* valuesField_ptr, 
									                        rgbaImage_t* imageField_ptr, 
															const COLOR::colorInterpolation_t* scheme_ptr);
	F_V2::texRetCode_st translateValuesToInterpolatedColors(const floats2Dfield_t* valuesField_ptr, 
									                        rgbaImage_t* imageField_ptr, 
															const COLOR::colorInterpolation_t* scheme_ptr);
	F_V2::texRetCode_st translateValuesToInterpolatedColors(const doubles2Dfield_t* valuesField_ptr, 
															rgbaImage_t* imageField_ptr, 
															const COLOR::colorInterpolation_t* scheme_ptr);

	typedef enum class imageType { PNG, JPG } imageType_t;

	//TODO: for now, "path" should include the slash. Change this and make it portable (possibly use fAux)
	//Quality affects quality only when type is JPG
	//Returns error code on errors or OK on success
	//If image_ptr is a field of floats or doubles, a temp greyscale image is used, using min and span to convert
	F_V2::imageFileRetCode_st saveImage(const generic2DfieldPtr_t* image_ptr, std::string filename, 
		                                imageType_t type, int quality = 100, std::string path = "",
		                                double min = 0, double span = 1);
}
