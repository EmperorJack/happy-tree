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

	// Build the particle system
	buildSystem();
}

FuzzyObject::~FuzzyObject() { }

void FuzzyObject::buildSystem() {
	for (int i = 0; i < maxParticles; i++) {
		particle p;
	  p.pos = vec3(rand() % 30 - 15, rand() % 30 - 15, rand() % 30 - 15);
	  particles.push_back(p);
	}
}

void FuzzyObject::updateSystem() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(rand() % 60 - 30, rand() % 60 - 30, rand() % 60 - 30) - particles[i].pos;
		particles[i].acc = clamp(normalize(particles[i].acc), -0.005f, 0.005f);

		particles[i].vel += particles[i].acc;
		particles[i].pos += particles[i].vel;
	}
}

void FuzzyObject::renderSystem() {
	// Set material properties
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, diffuse.dataPointer());
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular.dataPointer());
	glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	// Set drawing properties
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < particles.size(); i++) {
		particle p = particles[i];

		glPushMatrix();

		glTranslatef(p.pos.x, p.pos.y, p.pos.z);
		cgraSphere(pRadius, 4, 4);

		glPopMatrix();
	}
}
