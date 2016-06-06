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

FuzzyObject::FuzzyObject(Geometry *geometry) {
	g_geometry = geometry;

	setupDisplayList();

	//buildSystem(false);
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
	cgraSphere(p_radius, 8, 8);

	glEndList();
}

void FuzzyObject::buildSystemIncrement() {
	buildSystem(true);
}

void FuzzyObject::buildSystem(bool incremental) {
	// First pass
	while (!firstPassFinished && !stoppingCriteria()) {
		addParticle();
		updateSystem();

		if (incremental) return;
	}

	// Second pass
	int maxParticles = particles.size() - 1;
	while (particles.size() < maxParticles) {
		addParticle();
		updateSystem();

		if (incremental) return;
	}

	// Final pass
	while (!systemAtRest()) {
		updateSystem();

		if (incremental) return;
	}

	buildFinished = true;
}

bool FuzzyObject::stoppingCriteria() {
	if (particles.size() >= particleLimit) return true;

	if (collisionCount == particles.size() && particles.size() > 1) {

		for (int i = 0; i < stabilityUpdates; i++) {
			updateSystem();
		}

		if (collisionCount == particles.size()) {
			firstPassFinished = true;
			return true;
		}
	}

	return false;
}

// 	// Compute the total velocity of the system
// 	vec3 totalVelocity;
// 	for (int i = 0; i < particles.size(); i++) {
// 		totalVelocity += particles[i].vel;
// 	}

// 	bool stable = collisionCount == 0;

// 	// If the system regained stability last step but lost it again
// 	if (wasStable && !stable) {
// 		stabilityStepCount++;

// 		if (stabilityStepCount >= stabilityThreshold && length(totalVelocity) < p_velRange) {
// 			// System cannot maintain stability, it is overfull
// 			finished = true;
// 			return true;
// 		}
// 	} else {
// 		// Reset stability counter
// 		wasStable = false;
// 		stabilityStepCount = 0;
// 	}
// 	// TODO Make criteria for below when all particles are colliding
// 	// If the particles are repeatedly colliding
// 	if (!overfull && !stable && collisionCount >= lastCollisionCount) {
// 		// Move each particle slightly
// 		for (int i = 0; i < particles.size(); i++) {
// 			particles[i].pos.x += math::random(-manualShiftAmount, manualShiftAmount);
// 			particles[i].pos.y += math::random(-manualShiftAmount, manualShiftAmount);
// 			particles[i].pos.z += math::random(-manualShiftAmount, manualShiftAmount);
// 		}

// 		overfullStepCount++;
// 	} else {
// 		// Reset overfull counter
// 		overfullStepCount = 0;
// 	}

// 	// Check if the system is overfull
// 	if (!overfull && overfullStepCount >= overfullThreshold) {
// 		overfull = true;
// 		overfullStepCount = 0;
// 		cout << "overfull" << endl;
// 		// Remove a random particle
// 		particles.erase(particles.begin() + (int) math::random(0.0f, (float) particles.size()));
// 		//particles.pop_back();
// 	}

// 	// If the system was overfull check and it has regained stability
// 	if (overfull && stable) {
// 		// Enable particles to spawn again
// 		overfull = false;
// 		wasStable = true;
// 	}

// 	lastCollisionCount = collisionCount;

// 	return false;
// }

bool FuzzyObject::systemAtRest() {
	return true;
}

void FuzzyObject::addParticle() {
	// Do not add another particle if we reached the limit
	if (particles.size() >= particleLimit) return;

	particle p;
	p.pos = vec3(spawnPoint.x + math::random(-0.1f, 0.1f),
							 spawnPoint.y + math::random(-0.1f, 0.1f),
							 spawnPoint.z + math::random(-0.1f, 0.1f));

	p.acc = vec3(0.0f, 0.0f, 0.0f);

	// Random velocity generation
	p.vel = vec3(math::random(-1.0f, 1.0f) * p_velRange,
							 math::random(-1.0f, 1.0f) * p_velRange,
							 math::random(-1.0f, 1.0f) * p_velRange);

	particles.push_back(p);

	updateFacingTriangle(particles.size() - 1);
}

void FuzzyObject::updateSystem() {
	collisionCount = 0;

	// Reset each particles acceleration
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc = vec3(0.0f, 0.0f, 0.0f);
		particles[i].col = vec3(0.0f, 0.0f, 1.0f);
		particles[i].inCollision = false;

		// Check if the particle left the mesh
		float d = dot(particles[i].pos - particles[i].triangleIntersectionPos, -g_geometry->getSurfaceNormal(particles[i].triangleIndex));
		if (d < 0.0f || d == maxFloatVector.x) {
			particles[i].col = vec3(0.0f, 1.0f, 0.0f);
			particlesForDeletion.push_back(i);
		}
	}

	// Apply LJ physics based forces to the particle system
	applyParticleForces();

	// Apply forces to particles that collide with the geometry
	applyBoundaryForces();

	// Delete any particles marked for deletion
	for (int i = 0; i < particlesForDeletion.size(); i++) {
		particles.erase(particles.begin() + particlesForDeletion[i]);
	}
	particlesForDeletion.clear();

	// Update the particle positions and velocities
	for (int i = 0; i < particles.size(); i++) {
		particles[i].acc /= p_mass;
		particles[i].vel = clamp(particles[i].vel + particles[i].acc, -p_velRange, p_velRange);
		particles[i].pos += particles[i].vel;

		// If the particle accelerated it has potentially changed direction
		if (particles[i].acc.x != 0.0f && particles[i].acc.y != 0.0f && particles[i].acc.z != 0.0f) {

			// Recompute the particle facing triangle
			updateFacingTriangle(i);
		}

		if (particles[i].inCollision) collisionCount++;
	}
}

void FuzzyObject::applyParticleForces() {
	// For each pair of particles
	for (int i = 0; i < particles.size(); i++) {
		for (int j = i + 1; j < particles.size(); j++) {

			// If the particles are within the effect range of eachother we count this as a collision
			if (withinRange(particles[i].pos, particles[j].pos, e_effectRange)) {

				// Compute the distance between particles
				vec3 distVector = particles[i].pos - particles[j].pos;
				float dist = length(distVector);

				if (dist < 0.001f) continue; // Prevent dividing by 0 effects

				// Compute and apply the force both particles exert on eachother
				vec3 force = forceAtDistance(dist, distVector);
				particles[i].acc += force;
				particles[j].acc -= force;

				// Apply friction to the particle
				particles[i].vel *= particleCollisionFriction;
				particles[j].vel *= particleCollisionFriction;
			}
		}
	}
}

void FuzzyObject::applyBoundaryForces() {
	// For each particle
	for (int i = 0; i < particles.size(); i++) {

		// If the particle is colliding with the intersection point
		if (withinRange(particles[i].pos, particles[i].triangleIntersectionPos, p_boundaryRadius)) {
		 	// Bounce the particle off the triangle surface by reflecting it's velocity
			particles[i].vel = reflect(particles[i].vel, -(g_geometry->getSurfaceNormal(particles[i].triangleIndex))) * meshCollisionFriction;
			particles[i].acc = vec3(0.0f, 0.0f, 0.0f);

			// The particle is now facing the opposite direction so the facing triangle must be recomputed
			updateFacingTriangle(i);
		}
	}
}

vec3 FuzzyObject::forceAtDistance(float dist, vec3 distVector) {
	float a = (48 * e_strength / pow(e_lengthScale, 2));
	float b = pow(e_lengthScale / dist, 14);
	float c = 0.5f * pow(e_lengthScale / dist, 8);
	return a * (b - c) * distVector;
}

// Use the square distance to cut costs by avoiding square roots
bool FuzzyObject::withinRange(vec3 p1, vec3 p2, float range) {
	vec3 d = p1 - p2;
	return (d.x * d.x + d.y * d.y + d.z * d.z) < (range * range);
}

void FuzzyObject::updateFacingTriangle(int index) {
	vec3 newIntersectionPoint = vec3(maxFloatVector);
	float shortestLength = maxFloatVector.x;
	int triangleIndex = 0;

	// For each triangle
	for (int i = 0; i < g_geometry->triangleCount(); i++) {

		// Using the particle velocity as the direction vector
		// Compute the intersection point on the triangle
		vec3 intersectionPoint = (g_geometry->rayIntersectsTriangle(particles[index].pos, particles[index].vel, i));

		// Skip this triangle if no intersection occured
		if (intersectionPoint.x == maxFloatVector.x) continue;

		// If this is the closest intersection point yet
		if (withinRange(particles[index].pos, intersectionPoint, shortestLength)) {
			newIntersectionPoint = intersectionPoint;
			shortestLength = length(particles[index].pos - newIntersectionPoint);
			triangleIndex = i;
		}
	}

	// Assign the final closest intersection point
	particles[index].triangleIntersectionPos = newIntersectionPoint;
	particles[index].triangleIndex = triangleIndex;

	particles[index].col = vec3(1.0f, 0.0f, 0.0f);
	particles[index].inCollision = true;
}

void FuzzyObject::renderSystem() {
	// Draw the spawn position
	glDisable(GL_LIGHTING);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vec4(1.0f, 0.0f, 0.0f, 1.0f).dataPointer());
	glPointSize(4);

	glBegin(GL_POINTS);
	glVertex3f(spawnPoint.x, spawnPoint.y, spawnPoint.z);
	glEnd();

	glEnable(GL_LIGHTING);

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
		if (particleViewMode) {
			glCallList(p_displayList);
		} else {
			glBegin(GL_POINTS);
			glVertex3f(0, 0, 0);
			glEnd();

			// Draw the particle velocity
			glBegin(GL_LINES);
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(p.vel.x * 2, p.vel.y * 2, p.vel.z * 2);
			glEnd();
		}

		glPopMatrix();
	}
	glEnd();
}

int FuzzyObject::getParticleCount() {
	return particles.size();
}

void FuzzyObject::toggleParticleViewMode() {
	particleViewMode = !particleViewMode;
}
