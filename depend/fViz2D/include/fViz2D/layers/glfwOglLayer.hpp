#pragma once

#include "GLFW/glfw3.h" // Will drag system OpenGL headers

#include "utils/imageUtils.hpp"

GLFWwindow* initGlfwAndCreateWindow(GLFWerrorfun glfwErrorCallback, int width, int height, 
                                                     const char* title, bool vsync = true);
void shutDownGLFW(GLFWwindow* window);

void clearFrameBuffer(GLFWwindow* window, COLOR::rgbaF_t clearColor);