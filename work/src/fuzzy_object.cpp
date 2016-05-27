//---------------------------------------------------------------------------
// Particle Representation of a 3D Object
//
// By Jack Purvis
//
// Referenced paper:
// “Obtaining Fuzzy Representations of 3D Objects”, by Brent M. Dingle, November 2005.
// https://engineering.tamu.edu/media/697054/tamu-cs-tr-2005-11-6.pdf
//---------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "cgra_math.hpp"
#include "cgra_geometry.hpp"
#include "fuzzy_object.hpp"
#include "geometry.hpp"
#include "opengl.hpp"

using namespace std;
using namespace cgra;

FuzzyObject::FuzzyObject(Geometry *g) {
	geometry = g;

	setupDisplayList();

	cout << "Length scale: " << e_lengthScale << endl;
	cout << "Effect range: " << e_effectRange << endl;
}

FuzzyObject::~FuzzyObject() {}

// Setup the particle instance geometry
void FuzzyObject::setupDisplayList() {
	// Delete the old list if there is one
	if (p_displayList) glDeleteLists(p_displayList, 1);

	// Setup the new display list
	p_displayList = glGenLists(1);
	glNewList(p_displayList, GL_COMPILE);

	// Draw the geometry
	cgraSphere(p_radius, 3, 3);

	glEndList();
}

void FuzzyObject::buildSystem() {}

void FuzzyObject::buildIncremental() {
	// First pass
	if (!stoppingCriteria()) {
		addParticle();
		updateSystem();
		return;
	}

	// Second pass
	int maxParticles = particles.size() - 1;
	if (particles.size() < maxParticles) {
		addParticle();
		//updateSystem();
		return;
	}

	// Final pass
	if (!systemAtRest()) {
		updateSystem();
	}
}

void FuzzyObject::addParticle() {
	particle p;
	p.pos = spawnPoint;
	//p.pos = vec3(spawnPoint.x + math::random(-1.0f, 1.0f),
	//						 spawnPoint.y + math::random(-1.0f, 1.0f),
	//						 spawnPoint.z + math::random(-1.0f, 1.0f));

	p.acc = vec3(0.0f, 0.0f, 0.0f);

	// Random velocity generation
	p.vel = vec3(math::random(-1.0f, 1.0f) * p_velRange,
							 math::random(-1.0f, 1.0f) * p_velRange,
							 math::random(-1.0f, 1.0f) * p_velRange);

	p.col = vec3(math::random(0.0f, 1.0f),
							 math::random(0.0f, 1.0f),
							 math::random(0.0f, 1.0f));

	particles.push_back(p);
}

void FuzzyObject::updateSystem() {
	// For each particle
	for (int i = 0; i < particles.size(); i++) {

		vec3 forceVector = vec3(0.0f, 0.0f, 0.0f);

		// For each other particle
		for (int j = 0; j < particles.size(); j++) {

			if (i == j) continue;

			vec3 distVector = particles[i].pos - particles[j].pos;
			float dist = length(distVector);

			if (dist < e_effectRange && dist > e_lengthScale) {
				forceVector += forceAtDistance(dist, distVector);
			}
		}

		particles[i].acc += forceVector / p_mass;
	}

	// For each particle
	for (int i = 0; i < particles.size(); i++) {
		particles[i].vel += particles[i].acc;
		particles[i].pos += particles[i].vel;
	}
}

vec3 FuzzyObject::forceAtDistance(float dist, vec3 distVector) {
	float a = (48 * e_strength / pow(e_lengthScale, 2));
	float b = pow(e_lengthScale / dist, 14);
	float c = 0.5f * pow(e_lengthScale / dist, 8);
	return a * (b - c) * distVector;
}

bool FuzzyObject::stoppingCriteria() {
	return particles.size() >= particleLimit;
}

bool FuzzyObject::systemAtRest() {
	return false;
}

void FuzzyObject::renderSystem() {
	// Draw the spawn position
	glDisable(GL_LIGHTING);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vec4(1.0f, 0.0f, 0.0f, 1.0f).dataPointer());
	glPointSize(8);

	glBegin(GL_POINTS);
	glVertex3f(spawnPoint.x, spawnPoint.y, spawnPoint.z);
	glEnd();

	// Draw the effect radius
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vec4(0.0f, 0.0f, 1.0f, 1.0f).dataPointer());
	glLineWidth(6);

	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(e_effectRange, 0.0f, 0.0f);
	glEnd();

	glEnable(GL_LIGHTING);

	// Set particle material properties
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular.dataPointer());
	glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	// Set particle drawing properties
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
	for (int i = 0; i < particles.size(); i++) {
		particle p = particles[i];

		glPushMatrix();
		glTranslatef(p.pos.x, p.pos.y, p.pos.z);

		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, p.col.dataPointer());
		glCallList(p_displayList);
		//glBegin(GL_POINTS);
		//glVertex3f(p.pos.x, p.pos.y, p.pos.z);
		//glEnd();

		glPopMatrix();
	}
	glEnd();
}

int FuzzyObject::getParticleCount() {
	return particles.size();
}
