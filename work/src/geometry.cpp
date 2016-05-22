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

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "cgra_math.hpp"
#include "geometry.hpp"
#include "opengl.hpp"

using namespace std;
using namespace cgra;

Geometry::Geometry(string filename) {
	m_filename = filename;
	readOBJ(filename);
	if (m_triangles.size() > 0) {
		createDisplayListPoly();
		createDisplayListWire();
	}

	// Default material setting
	m_material.ambient = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	m_material.diffuse = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	m_material.specular = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	m_material.shininess = 0.0f;
	m_material.emission = vec4(0.0f, 0.0f, 0.0f, 0.0f);
}

Geometry::~Geometry() { }

void Geometry::readOBJ(string filename) {

	// Make sure our geometry information is cleared
	m_points.clear();
	m_uvs.clear();
	m_normals.clear();
	m_triangles.clear();

	// Load dummy points because OBJ indexing starts at 1 not 0
	m_points.push_back(vec3(0,0,0));
	m_uvs.push_back(vec2(0,0));
	m_normals.push_back(vec3(0,0,1));

	ifstream objFile(filename);

	if(!objFile.is_open()) {
		cerr << "Error reading " << filename << endl;
		throw runtime_error("Error :: could not open file.");
	}

	cout << "Reading file " << filename << endl;

	// good() means that failbit, badbit and eofbit are all not set
	while(objFile.good()) {

		// Pull out line from file
		string line;
		std::getline(objFile, line);
		istringstream objLine(line);

		// Pull out mode from line
		string mode;
		objLine >> mode;

		// Reading like this means whitespace at the start of the line is fine
		// attempting to read from an empty string/line will set the failbit
		if (!objLine.fail()) {

			if (mode == "v") {
				vec3 v;
				objLine >> v.x >> v.y >> v.z;
				m_points.push_back(v);

			} else if(mode == "vn") {
				vec3 vn;
				objLine >> vn.x >> vn.y >> vn.z;
				m_normals.push_back(vn);

			} else if(mode == "vt") {
				vec2 vt;
				objLine >> vt.x >> vt.y;
				m_uvs.push_back(vt);

			} else if(mode == "f") {

				vector<vertex> verts;
				while (objLine.good()){
					vertex v;

					objLine >> v.p;		// Scan in position index

					// Check if the obj file contains normal data
					if (m_normals.size() > 1) {
						// obj face is formatted as v/vt/vn or v//vn
						objLine.ignore(1);	// Ignore the '/' character

						if (m_uvs.size() > 1) { // Check if the obj file contains uv data
							objLine >> v.t;		// Scan in uv (texture coord) index
						}
						objLine.ignore(1);	// Ignore the '/' character

						objLine >> v.n;		// Scan in normal index
					} // Else the obj face is formatted just as v

					verts.push_back(v);
				}

				// If we have 3 verticies, construct a triangle
				if(verts.size() >= 3){
					triangle tri;
					tri.v[0] = verts[0];
					tri.v[1] = verts[1];
					tri.v[2] = verts[2];
					m_triangles.push_back(tri);
				}
			}
		}
	}

	cout << "Reading OBJ file is DONE." << endl;
	cout << m_points.size()-1 << " points" << endl;
	cout << m_uvs.size()-1 << " uv coords" << endl;
	cout << m_normals.size()-1 << " normals" << endl;
	cout << m_triangles.size() << " faces" << endl;

	// If we didn't have any normals, create them
	if (m_normals.size() <= 1) createNormals();
}

void Geometry::createNormals() {
	// Maintain a collection of each vertex normal
	vector<vec3> vertex_normals(m_points.size());

	for (int i = 0; i < m_triangles.size(); i++) {
		vec3 v1 = m_points[m_triangles[i].v[0].p];
		vec3 v2 = m_points[m_triangles[i].v[1].p];
		vec3 v3 = m_points[m_triangles[i].v[2].p];

		// Compute the surface normal from the triangle vertices
		vec3 n = cgra::cross(v2 - v1, v3 - v1);

		// Add the surface normal to each associated vertex normal
		for (int j = 0; j < 3; j++) {
			vertex_normals[m_triangles[i].v[j].p] += n;
		}
	}

	for (int i = 0; i < m_triangles.size(); i++) {
		for (int j = 0; j < 3; j++) {
			// Push the normalized vertex normal to the model normals
			vec3 v_normal = vertex_normals[m_triangles[i].v[j].p];
			m_normals.push_back(cgra::normalize(v_normal));

			// Associate the current vertex with the new vertex normal
			m_triangles[i].v[j].n = m_normals.size() - 1;
		}
	}
}

void Geometry::createDisplayListPoly() {
	// Delete old list if there is one
	if (m_displayListPoly) glDeleteLists(m_displayListPoly, 1);

	// Create a new list
	cout << "Creating Poly Geometry" << endl;
	m_displayListPoly = glGenLists(1);
	glNewList(m_displayListPoly, GL_COMPILE);

	displayTriangles();

	glEndList();
	cout << "Finished creating Poly Geometry" << endl;
}

void Geometry::createDisplayListWire() {
	// Delete old list if there is one
	if (m_displayListWire) glDeleteLists(m_displayListWire, 1);

	// Create a new list
	cout << "Creating Wire Geometry" << endl;
	m_displayListWire = glGenLists(1);
	glNewList(m_displayListWire, GL_COMPILE);

	displayTriangles();

	glEndList();
	cout << "Finished creating Wire Geometry" << endl;
}

void Geometry::displayTriangles() {
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < m_triangles.size(); i++) {
		for (int j = 0; j < 3; j++) {
			vec3 p = m_points[m_triangles[i].v[j].p];
			vec2 t = m_uvs[m_triangles[i].v[j].t];
			vec3 n = m_normals[m_triangles[i].v[j].n];

			glNormal3f(n.x, n.y, n.z);
			glTexCoord2f(t.x * 4.0f, t.y * 4.0f);
			glVertex3f(p.x, p.y, p.z);
		}
	}
	glEnd();
}

void Geometry::setPosition(vec3 position) {
	m_position = position;
}

void Geometry::setMaterial(vec4 ambient, vec4 diffuse, vec4 specular, float shininess, vec4 emission) {
	m_material.ambient = ambient;
	m_material.diffuse = diffuse;
	m_material.specular = specular;
	m_material.shininess = shininess;
	m_material.emission = emission;
}

void Geometry::renderGeometry() {
	glPushMatrix();

	// Translate to the object position
	glTranslatef(m_position.x, m_position.y, m_position.z);

	// Apply the object material properties
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m_material.diffuse.dataPointer());
	glMaterialfv(GL_FRONT, GL_SPECULAR, m_material.specular.dataPointer());
	glMaterialfv(GL_FRONT, GL_SHININESS, &m_material.shininess);
	glMaterialfv(GL_FRONT, GL_EMISSION, m_material.emission.dataPointer());

	glShadeModel(GL_SMOOTH);

	if (wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glCallList(m_displayListWire);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glCallList(m_displayListPoly);
	}

	glPopMatrix();
}

void Geometry::toggleWireframe() {
	wireframe = !wireframe;
}
