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

#include "cgra_math.hpp"
#include "opengl.hpp"

#include <stb_image.h>


class Image {
public:
	int w, h, n;
	std::vector<unsigned char> data;

	Image(int w_, int h_, int n_) : w(w_), h(h_), n(n_), data(w*h*n, 0) {}

	Image(const std::string &filepath) {
		unsigned char *stbi_data = stbi_load(filepath.c_str(), &w, &h, &n, 0);
		if (stbi_data == NULL) throw std::runtime_error("Error: Failed to load image " + filepath + " : file doesn't exist or is an unsupported format.");
		if (n > 4) throw std::runtime_error("Error: Failed to load image " + filepath + " : greater than 4 channels not supported.");
		data.assign(stbi_data, stbi_data + (w*h*n));
		stbi_image_free(stbi_data);
	}

	Image(const Image &) = default;
	Image & operator=(const Image &) = default;
	Image(Image &&) = default;
	Image & operator=(Image &&) = default;


	// Use to get the appropriate GL format for data
	GLenum glFormat() const {
		switch (n) {
		case 1: return GL_R;
		case 2: return GL_RG;
		case 3: return GL_RGB;
		case 4: return GL_RGBA;
		default: return GL_RGB; // TODO list
		}
	}

	// Use to get a GL friendly pointer to the data
	unsigned char * dataPointer() { return &data[0]; }
	const unsigned char * dataPointer() const { return &data[0]; }

	Image subsection(int xoffset, int yoffset, int width, int height) {
		Image r(width, height, n);

		for (int y = 0; y < height; y++) {
			if ((y + yoffset) >= h) continue;
			for (int x = 0; x < width; x++) {
				if ((x + xoffset) >= w) continue;
				for (int i = 0; i < n; i++) {
					r.data[(y*width*n) + (x*n) + i] =
						data[((y + yoffset)*w*n) + ((x + xoffset)*n) + i];
				}
			}
		}
		return r;
	}
};