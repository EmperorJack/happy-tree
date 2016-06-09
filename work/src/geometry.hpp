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
		Geometry(std::string, float texScale = 1.0f);
		Geometry(std::vector<cgra::vec3> points,
				   std::vector<cgra::vec3> normals,
				   std::vector<cgra::vec2> uvs,
				   std::vector<triangle> triangles,
				   bool);
		~Geometry();

		cgra::vec3 getPosition();
		void setPosition(cgra::vec3);
		void setMaterial(cgra::vec4, cgra::vec4, cgra::vec4, float, cgra::vec4);
		cgra::vec3 rayIntersectsTriangle(cgra::vec3, cgra::vec3, int);
		bool pointInsideMesh(cgra::vec3);
		void renderGeometry(bool);
		int triangleCount();
		cgra::vec3 getSurfaceNormal(int);
		cgra::vec3 getOrigin();

	private:
		std::string m_filename;

		// Fields for storing raw obj information
		std::vector<cgra::vec3> m_points;
		std::vector<cgra::vec2> m_uvs;
		std::vector<cgra::vec3> m_normals;
		std::vector<triangle> m_triangles;
		std::vector<cgra::vec3> m_surfaceNormals;

		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);
		material m_material;
		float m_textureScale = 1.0f;

		// The vector to return if no ray intersection is found
		cgra::vec3 noIntersectionVector = cgra::vec3(std::numeric_limits<float>::max(), 0.0f, 0.0f);

		// IDs for the display list to render
		GLuint m_displayListPoly = 0;
		GLuint m_displayListWire = 0;

		void readOBJ(std::string);
		void createNormals();
		void createSurfaceNormals();
		void createDisplayListPoly();
		void createDisplayListWire();
		void displayTriangles();
};
