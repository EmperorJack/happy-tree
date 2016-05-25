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

	cout << effectRange << endl;
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
		updateSystem();
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
	p.acc = vec3(0.0f, 0.0f, 0.0f);

	// Random velocity generation
	p.vel = vec3(math::random(-1.0f, 1.0f) * p_velRange, math::random(-1.0f, 1.0f) * p_velRange, math::random(-1.0f, 1.0f) * p_velRange);

	particles.push_back(p);
}

void FuzzyObject::updateSystem() {
	// For each particle
	for (int i = 0; i < particles.size(); i++) {

		// For each other particle
		for (int j = i + 1; j < particles.size(); j++) {

			float dist = distance(particles[i].pos, particles[j].pos);
			if (dist < effectRange) {
				//
			} else {
				particles[i].acc = vec3(0.0f, 0.0f, 0.0f);
			}
		}

		particles[i].vel += particles[i].acc;
		particles[i].pos += particles[i].vel;
	}
}

bool FuzzyObject::stoppingCriteria() {
	return false;
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

	glEnable(GL_LIGHTING);

	// Set particle material properties
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, diffuse.dataPointer());
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular.dataPointer());
	glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	// Set particle drawing properties
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < particles.size(); i++) {
		particle p = particles[i];

		glPushMatrix();

		glTranslatef(p.pos.x, p.pos.y, p.pos.z);
		glCallList(p_displayList);

		glPopMatrix();
	}
}

int FuzzyObject::getParticleCount() {
	return particles.size();
}
