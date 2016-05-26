//---------------------------------------------------------------------------
//
// Copyright (c) 2016 Taehyun Rhee, Joshua Scott, Ben Allen
//
// This software is provided 'as-is' for assignment of COMP308 in ECS,
// Victoria University of Wellington, without any express or implied warranty. 
// In no event will the authors be held liable for any damages arising from
// the use of this software.
//
// The contents of this file may not be copied or duplicated in any form
// without the prior permission of its owner.
//
//----------------------------------------------------------------------------

#pragma once

#include "imgui.h"
#include "opengl.hpp"

namespace cgra {
	namespace SimpleGUI {
		void mouseButtonCallback(GLFWwindow*, int button, int action, int /*mods*/);
		void scrollCallback(GLFWwindow*, double /*xoffset*/, double yoffset);
		void keyCallback(GLFWwindow*, int key, int /*scancode*/, int action, int mods);
		void charCallback(GLFWwindow*, unsigned int c);
		bool init(GLFWwindow* window, bool install_callbacks=false);
		void shutdown();
		void newFrame();
		void render();
	}
}