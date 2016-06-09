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
	// Create a particle system out of given points
	for (int i = 0; i < points.size(); i++) {
		particle p;
		p.pos = points[i];
		p.vel = vec3(0.0f, 0.0f, 0.0f);
		p.acc = vec3(0.0f, 0.0f, 0.0f);
		particles.push_back(p);
	}

	setupDisplayList();
}

void ParticleSystem::update() {
	// For each particle
	for (int i = 0; i < particles.size(); i++) {
		particles[i].vel = particles[i].vel + particles[i].acc;
		particles[i].pos += particles[i].vel;

		vec3 actualPos = particles[i].pos;
		if (actualPos.y - p_radius < 0.0f) {
			particles[i].vel.y *= -1;
			particles[i].vel *= 0.9f;
			particles[i].pos.y += p_radius;
		}
	}
}

void ParticleSystem::render() {
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

		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, p.col.dataPointer());

		// Draw the particle
		glCallList(p_displayList);

		glPopMatrix();
	}

	glPopMatrix();
}

// Setup the particle instance geometry
void ParticleSystem::setupDisplayList() {
	// Delete the old list if there is one
	if (p_displayList) glDeleteLists(p_displayList, 1);

	// Setup the new display list
	p_displayList = glGenLists(1);
	glNewList(p_displayList, GL_COMPILE);

	// Draw the geometry
	cgraSphere(p_radius, 8, 8);

	glEndList();
}

void ParticleSystem::drop() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(0.0f, -0.00981f, 0.0f);
		particles[i].vel = vec3(math::random(-1.0f, 1.0f) * p_velRange / 2.0f,
			                      math::random(-0.01f, 0.0f),
			                      math::random(-1.0f, 1.0f) * p_velRange / 2.0f);
	}
}

void ParticleSystem::explode() {
	for (int i = 0; i < particles.size(); i++) {
		particles[i].vel = vec3(math::random(-1.0f, 1.0f) * p_velRange * 10.0f,
			                      math::random(-1.0f, 1.0f) * p_velRange * 10.0f,
			                      math::random(-1.0f, 1.0f) * p_velRange * 10.0f);
	}
}
