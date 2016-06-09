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

Geometry::Geometry(vector<vec3> points, vector<vec3> normals, vector<vec2> uvs, vector<triangle> triangles, bool genSurfaceNormals) {
	m_points = points;
	m_normals = normals;
	m_triangles = triangles;
	m_uvs = uvs;

	// Load a dummy point for the UV's
	//m_uvs.push_back(vec2(0,0));

	// Create the surface normals for every triangle
	if (genSurfaceNormals) createSurfaceNormals();

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

	// Create the surface normals for every triangle
	createSurfaceNormals();

	// If we didn't have any vertex normals, create them
	if (m_normals.size() <= 1) createNormals();

	cout << "Reading OBJ file is DONE." << endl;
	cout << m_points.size()-1 << " points" << endl;
	cout << m_uvs.size()-1 << " uv coords" << endl;
	cout << m_normals.size()-1 << " normals" << endl;
	cout << m_triangles.size() << " faces" << endl;
}

void Geometry::createNormals() {
	// Maintain a collection of each vertex normal
	vector<vec3> vertex_normals(m_points.size());

	for (int i = 0; i < m_triangles.size(); i++) {
		vec3 v1 = m_points[m_triangles[i].v[0].p];
		vec3 v2 = m_points[m_triangles[i].v[1].p];
		vec3 v3 = m_points[m_triangles[i].v[2].p];

		// Compute the surface normal from the triangle vertices
		vec3 n = cross(v2 - v1, v3 - v1);

		// Add the surface normal to each associated vertex normal
		for (int j = 0; j < 3; j++) {
			vertex_normals[m_triangles[i].v[j].p] += n;
		}
	}

	for (int i = 0; i < m_triangles.size(); i++) {
		for (int j = 0; j < 3; j++) {
			// Push the normalized vertex normal to the model normals
			vec3 v_normal = vertex_normals[m_triangles[i].v[j].p];
			m_normals.push_back(normalize(v_normal));

			// Associate the current vertex with the new vertex normal
			m_triangles[i].v[j].n = m_normals.size() - 1;
		}
	}
}

void Geometry::createSurfaceNormals() {
	// Compute the surface normal of each triangle
	for (int i = 0; i < m_triangles.size(); i++) {
		vec3 surfaceNormal = cross(m_points[m_triangles[i].v[1].p] - m_points[m_triangles[i].v[0].p],
			 																	m_points[m_triangles[i].v[2].p] - m_points[m_triangles[i].v[0].p]);
		m_surfaceNormals.push_back(normalize(surfaceNormal));
	}
}

void Geometry::createDisplayListPoly() {
	// Delete old list if there is one
	if (m_displayListPoly) glDeleteLists(m_displayListPoly, 1);

	// Create a new list
	//cout << "Creating Poly Geometry" << endl;
	m_displayListPoly = glGenLists(1);
	glNewList(m_displayListPoly, GL_COMPILE);

	displayTriangles();

	glEndList();
	//cout << "Finished creating Poly Geometry" << endl;
}

void Geometry::createDisplayListWire() {
	// Delete old list if there is one
	if (m_displayListWire) glDeleteLists(m_displayListWire, 1);

	// Create a new list
	//cout << "Creating Wire Geometry" << endl;
	m_displayListWire = glGenLists(1);
	glNewList(m_displayListWire, GL_COMPILE);

	displayTriangles();

	glEndList();
	//cout << "Finished creating Wire Geometry" << endl;
}

void Geometry::displayTriangles() {
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < m_triangles.size(); i++) {
		for (int j = 0; j < 3; j++) {
			vec3 p = m_points[m_triangles[i].v[j].p];
			vec2 t = m_uvs[m_triangles[i].v[j].t];
			vec3 n = m_normals[m_triangles[i].v[j].n];

			glNormal3f(n.x, n.y, n.z);
			glTexCoord2f(t.x, t.y);
			glVertex3f(p.x, p.y, p.z);
		}
	}
	glEnd();
}

vec3 Geometry::getPosition() {
	return m_position;
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

vec3 Geometry::rayIntersectsTriangle(vec3 p, vec3 d, int triIndex) {
	// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

	triangle tri = m_triangles[triIndex];
	vec3 v0 = m_points[tri.v[0].p];
	vec3 v1 = m_points[tri.v[1].p];
	vec3 v2 = m_points[tri.v[2].p];

	vec3 e1, e2, h, s, q;
	float a,f,u,v;
	e1 = v1 - v0;
	e2 = v2 - v0;

	h = cross(d, e2);
	a = dot(e1, h);

	if (a > -0.00001f && a < 0.00001f) return noIntersectionVector;

	f = 1/a;
	s = p - v0;
	u = f * (dot(s,h));

	if (u < 0.0f || u > 1.0f) return noIntersectionVector;

	q = cross(s, e1);
	v = f * dot(d,q);

	if (v < 0.0f || u + v > 1.0f) return noIntersectionVector;

	// Compute t to find out where the intersection point is on the line
	float t = f * dot(e2,q);

	// Ray intersection occured
	if (t > 0.00001f) {
		return vec3(((1 - u - v) * v0) + (u * v1) + (v * v2));
	}

	// Line intersection occured but not a ray intersection
	return noIntersectionVector;
}

bool Geometry::pointInsideMesh(vec3 point) {
	int intersectionCount = 0;

	vec3 direction = vec3(0, 0, 1);

	for (int i = 0; i < m_triangles.size(); i++) {
		if (rayIntersectsTriangle(point, direction, i).x != noIntersectionVector.x) {
			intersectionCount++;
		}
	}

	// An odd number of intersections means the point is inside the mesh
	return intersectionCount % 2;
}

void Geometry::renderGeometry(bool wireframe) {
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
		glLineWidth(1);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glCallList(m_displayListWire);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glCallList(m_displayListPoly);
	}

	//Debug code for drawing the surface normals
	// for (int i = 0; i < m_surfaceNormals.size(); i++) {
	// 	glPushMatrix();
	// 	vec3 triPos = (m_points[m_triangles[i].v[0].p] + m_points[m_triangles[i].v[1].p] + m_points[m_triangles[i].v[2].p]) / 3.0f;
	// 	glTranslatef(triPos.x, triPos.y, triPos.z);

	// 	glBegin(GL_LINES);
	// 	glVertex3f(0.0f, 0.0f, 0.0f);
	// 	glVertex3f(m_surfaceNormals[i].x * 0.5f, m_surfaceNormals[i].y * 0.5f, m_surfaceNormals[i].z * 0.5f);
	// 	glEnd();

	// 	glBegin(GL_POINTS);
	// 	glVertex3f(m_surfaceNormals[i].x * 0.5f, m_surfaceNormals[i].y * 0.5f, m_surfaceNormals[i].z * 0.5f);
	// 	glEnd();

	// 	glPopMatrix();
	// }

	glPopMatrix();
}

int Geometry::triangleCount() {
  return m_triangles.size();
}

vec3 Geometry::getSurfaceNormal(int index) {
	return m_surfaceNormals[index];
}

vec3 Geometry::getOrigin() {
	vec3 average = vec3(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < m_points.size(); i++) {
		average += m_points[i];
	}
	return average / m_points.size();
}
