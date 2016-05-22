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
	  //p.vel = vec3((rand() % 100 / 100.0f - 0.5f) / 5.0f, (rand() % 100 / 100.0f - 0.5f) / 5.0f, (rand() % 100 / 100.0f - 0.5f) / 5.0f);
		particles.push_back(p);
	}
}

void FuzzyObject::updateSystem() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(rand() % 30 - 15, rand() % 30 - 15, rand() % 30 - 15) - particles[i].pos;
		particles[i].acc = clamp(normalize(particles[i].acc), -0.01f, 0.01f);

		particles[i].vel += particles[i].acc;
		particles[i].pos += particles[i].vel;
	}
}

void FuzzyObject::renderSystem() {
	for (int i = 0; i < particles.size(); i++) {
		particle p = particles[i];

		glPushMatrix();

		glTranslatef(p.pos.x, p.pos.y, p.pos.z);
		cgraSphere(pRadius, 8, 8);

		glPopMatrix();
	}
}
