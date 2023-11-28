#pragma once

//TODO: review: some stuff here probably could be made more general to support more types (template?)

#include "GLFW/glfw3.h"

#include "fViz2D/utils/imageUtils.hpp"
#include "fViz2D/API/returnCodes.hpp"

namespace TEX {

	//Holds data to identify and modify a texture on the GPU (using openGL 3.3+)
	typedef struct textureID_t {
		GLuint ID = 0;
		bool initialized = false;
		int width, height;

		bool createOGLtexID(GLenum minFilter = GL_LINEAR, GLenum magFilter = GL_LINEAR);
		~textureID_t() { glDeleteTextures(1, &ID); }
	} textureID_t;

	//Creates an OpenGL texture identifier and associates it with the rgbaTextureID_t
	bool createOGLtex(GLenum minFilter, GLenum magFilter, GLuint* ID_ptr);
	
	///Texture Loading:

	//Loads image file into an OpenGL texture with common settings. Returns an unitialized textureID_t on failure
	TEX::textureID_t load4channelTextureFromFile(const char* filename);	

	//Loads a generic data field into a texture. Returns BAD_DYNAMIC_DATA_FORMAT if format is not supported.
	//See specific load functions bellow for the supported formats and their behavior.
	F_V2::rendererRetCode_st loadTextureFromGeneric2DfieldPtr(IMG::generic2DfieldPtr_t* dataPtr_ptr, 
		                                                              TEX::textureID_t* texture_ptr);

	//Loads an IMG::rgbaImage_t image into a texture to be rendered.
	//If the textureID was unitialized, it's initialized first, creating an OpenGL texture identifier.
	void load4channelTextureFromRgbaImage(IMG::rgbaImage_t* image_ptr, TEX::textureID_t* texture_ptr);

	//Loads IMG::floats2Dfield_t field into a red-channel only texture of floats to be rendered.
	//If the textureID was unitialized, it's initialized first, creating an OpenGL texture identifier.
	void loadR32FtextureFromFloats(IMG::floats2Dfield_t* field_ptr, textureID_t* texture_ptr);
}