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

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "opengl.hpp"

struct vertex {
	int p = 0; // index for point in m_points
	int t = 0; // index for uv in m_uvs
	int n = 0; // index for normal in m_normals
};

struct triangle {
	vertex v[3]; //requires 3 verticies
};

struct material {
	cgra::vec4 ambient;
	cgra::vec4 diffuse;
	cgra::vec4 specular;
	float shininess;
	cgra::vec4 emission;
};

class Geometry {

	public:
		Geometry(std::string);
		~Geometry();

		void setPosition(cgra::vec3);
		void setMaterial(cgra::vec4, cgra::vec4, cgra::vec4, float, cgra::vec4);
		bool pointInsideMesh(cgra::vec3);
		void renderGeometry();
		void toggleWireframe();

	private:
		std::string m_filename;

		// Fields for storing raw obj information
		std::vector<cgra::vec3> m_points;
		std::vector<cgra::vec2> m_uvs;
		std::vector<cgra::vec3> m_normals;
		std::vector<triangle> m_triangles;

		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);
		material m_material;
		bool wireframe = false;

		// IDs for the display list to render
		GLuint m_displayListPoly = 0;
		GLuint m_displayListWire = 0;

		void readOBJ(std::string);
		void createNormals();
		void createDisplayListPoly();
		void createDisplayListWire();
		void displayTriangles();
		bool rayIntersectsTriangle(cgra::vec3, cgra::vec3, cgra::vec3, cgra::vec3, cgra::vec3);
};
