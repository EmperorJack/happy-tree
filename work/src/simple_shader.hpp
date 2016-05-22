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

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "opengl.hpp"

namespace cgra {


	class shader_error : public std::runtime_error {
	public:
		explicit shader_error(const std::string &what_ = "Generic shader error.") : std::runtime_error(what_) { }
	};

	class shader_type_error : public shader_error {
	public:
		explicit shader_type_error(const std::string &what_ = "Bad shader type.") : shader_error(what_) { }
	};

	class shader_compile_error : public shader_error {
	public:
		explicit shader_compile_error(const std::string &what_ = "Shader compilation failed.") : shader_error(what_) { }
	};

	class shader_link_error : public shader_error {
	public:
		explicit shader_link_error(const std::string &what_ = "Shader program linking failed.") : shader_error(what_) { }
	};

	inline void printShaderInfoLog(GLuint obj) {
		int infologLength = 0;
		int charsWritten = 0;
		glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infologLength);
		if (infologLength > 1) {
			std::vector<char> infoLog(infologLength);
			glGetShaderInfoLog(obj, infologLength, &charsWritten, &infoLog[0]);
			std::cout << "SimpleShader : " << "SHADER :\n" << &infoLog[0] << std::endl;
		}
	}

	inline void printProgramInfoLog(GLuint obj) {
		int infologLength = 0;
		int charsWritten = 0;
		glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);
		if (infologLength > 1) {
			std::vector<char> infoLog(infologLength);
			glGetProgramInfoLog(obj, infologLength, &charsWritten, &infoLog[0]);
			std::cout << "SimpleShader : " << "PROGRAM :\n" << &infoLog[0] << std::endl;
		}
	}

	inline GLuint compileShader(GLenum type, const std::string &text) {
		GLuint shader = glCreateShader(type);
		const char *text_c = text.c_str();
		glShaderSource(shader, 1, &text_c, nullptr);
		glCompileShader(shader);
		GLint compile_status;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);
		if (!compile_status) {
			printShaderInfoLog(shader);
			throw shader_compile_error();
		}
		// always print, so we can see warnings
		printShaderInfoLog(shader);
		return shader;
	}

	inline void linkShaderProgram(GLuint prog) {
		glLinkProgram(prog);
		GLint link_status;
		glGetProgramiv(prog, GL_LINK_STATUS, &link_status);
		if (!link_status) {
			printProgramInfoLog(prog);
			throw shader_link_error();
		}
		// always print, so we can see warnings
		printProgramInfoLog(prog);
	}

	inline GLuint makeShaderProgram(const std::vector<GLenum> &stypes, const std::vector<std::string> &sources) {
		if (stypes.size() != sources.size()) {
			throw std::runtime_error("Error: stypes and shader sources, vector size mismatch");
		}

		GLuint prog = glCreateProgram();

		for (size_t i = 0; i < stypes.size(); ++i) {
			auto shader = compileShader(stypes[i], sources[i]);
			glAttachShader(prog, shader);
		}

		linkShaderProgram(prog);
		std::cout << "SimpleShader : " << "Shader program compiled and linked successfully" << std::endl;
		return prog;
	}

	inline GLuint makeShaderProgramFromFile(const std::vector<GLenum> &stypes, const std::vector<std::string> &sourcefiles) {
		std::vector<std::string> sources;
		for (std::string filename : sourcefiles) {
			std::ifstream fileStream(filename);

			if (!fileStream) {
				throw std::runtime_error("Error: Could not locate and open file " + filename);
			}

			std::stringstream buffer;
			buffer << fileStream.rdbuf();
			sources.push_back(buffer.str());
		}

		return makeShaderProgram(stypes, sources);
	}

	inline GLuint makeShaderProgram(const std::string &profile, const std::vector<GLenum> &stypes, const std::string &source) {
		GLuint prog = glCreateProgram();

		auto get_define = [](GLenum stype) {
			switch (stype) {
			case GL_VERTEX_SHADER:
				return "_VERTEX_";
			case GL_GEOMETRY_SHADER:
				return "_GEOMETRY_";
			case GL_TESS_CONTROL_SHADER:
				return "_TESS_CONTROL_";
			case GL_TESS_EVALUATION_SHADER:
				return "_TESS_EVALUATION_";
			case GL_FRAGMENT_SHADER:
				return "_FRAGMENT_";
			default:
				return "_INVALID_SHADER_TYPE_";
			}
		};

		for (auto stype : stypes) {
			std::ostringstream oss;
			oss << "#version " << profile << std::endl;
			oss << "#define " << get_define(stype) << std::endl;
			oss << source;
			auto shader = compileShader(stype, oss.str());
			glAttachShader(prog, shader);
		}

		linkShaderProgram(prog);
		std::cout << "SimpleShader : " << "Shader program compiled and linked successfully" << std::endl;
		return prog;
	}

	inline GLuint makeShaderProgramFromFile(const std::string &profile, const std::vector<GLenum> &stypes, const std::string &sourcefile) {
		std::ifstream fileStream(sourcefile);

		if (!fileStream) {
			throw std::runtime_error("Error: Could not locate and open file " + sourcefile);
		}

		std::stringstream buffer;
		buffer << fileStream.rdbuf();

		return makeShaderProgram(profile, stypes, buffer.str());
	}	
}