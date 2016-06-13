//---------------------------------------------------------------------------
// General Particle System
//
// By Jack Purvis
//---------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "cgra_math.hpp"
#include "cgra_geometry.hpp"
#include "opengl.hpp"
#include "particle_system.hpp"

using namespace std;
using namespace cgra;

ParticleSystem::ParticleSystem(vector<vec3> points) {
	// Create a particle system out of the given points
	for (int i = 0; i < points.size(); i++) {
		particle p;
		p.original_pos = vec3(points[i]);
		p.pos = points[i];
		p.vel = vec3(0.0f, 0.0f, 0.0f);
		p.acc = vec3(0.0f, 0.0f, 0.0f);
		p.col = vec3(1.0f, 1.0f, 1.0f);
		particles.push_back(p);
	}

	setupDisplayList();
}

ParticleSystem::~ParticleSystem() {}

void ParticleSystem::update() {
	// For each particle
	for (int i = 0; i < particles.size(); i++) {

		// Update the particle velocity and position
		particles[i].vel = clamp(particles[i].vel + particles[i].acc, -p_maxVel, p_maxVel);
		particles[i].pos += particles[i].vel;

		// Perform simple bouncing off of the floor plane
		vec3 actualPos = particles[i].pos;
		if (actualPos.y - p_radius < 0.0f) {
			particles[i].vel.y *= -1;
			particles[i].vel *= 0.9f;
			particles[i].pos.y += p_radius;
		}
	}

	// Interpolate the current colour based on the animation state
	float lerp = animationStep / float(animationLength);
	currentColour = mix(startColour, endColour, lerp);
	animationStep = min(animationStep + 1, animationLength);
}

void ParticleSystem::render() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPushMatrix();

	// Set particle material properties
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular.dataPointer());
	glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	// Set particle drawing properties
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glLineWidth(1);

	for (int i = 0; i < particles.size(); i++) {
		particle p = particles[i];

		glPushMatrix();
		glTranslatef(p.pos.x, p.pos.y, p.pos.z);

		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColour.dataPointer());

		// Draw the particle
		glCallList(p_displayList);

		glPopMatrix();
	}

	glPopMatrix();

	glDisable(GL_BLEND);
}

// Setup the particle instance geometry
void ParticleSystem::setupDisplayList() {
	// Delete the old list if there is one
	if (p_displayList) glDeleteLists(p_displayList, 1);

	// Setup the new display list
	p_displayList = glGenLists(1);
	glNewList(p_displayList, GL_COMPILE);

	// Draw the geometry
	cgraSphere(p_radius, 6, 6);

	glEndList();
}

// Make each particle fall from it's current position
void ParticleSystem::drop() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(0.0f, -0.00981f, 0.0f);
		particles[i].vel = vec3(math::random(-1.0f, 1.0f) * p_velRange / 2.0f,
			                      math::random(-0.01f, 0.0f),
			                      math::random(-1.0f, 1.0f) * p_velRange / 2.0f);
	}
}

// Explode each particle in all directions
void ParticleSystem::explode() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(0.0f, -0.00981f, 0.0f);
		particles[i].vel = vec3(math::random(-1.0f, 1.0f) * p_velRange * 10.0f,
			                      math::random(-1.0f, 1.0f) * p_velRange * 10.0f,
			                      math::random(-1.0f, 1.0f) * p_velRange * 10.0f);
	}
}

// Attempt to blow the particle system away in the direction of a given wind vector
void ParticleSystem::blowAway(vec3 direction) {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(0.0f, -0.000981f, 0.0f);
		particles[i].acc += direction * math::random(1.0f, 5.0f);
		particles[i].vel = direction;
	}
}

// Reset the particle system to it's initial state
void ParticleSystem::resetParticles() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].pos = particles[i].original_pos;
		particles[i].vel = vec3(0.0f, 0.0f, 0.0f);
		particles[i].acc = vec3(0.0f, 0.0f, 0.0f);
	}
}
