//---------------------------------------------------------------------------
// Particle Representation of a 3D Object
//
// By Jack Purvis
//
// Referenced paper:
// “Obtaining Fuzzy Representations of 3D Objects”, by Brent M. Dingle, November 2005.
// https://engineering.tamu.edu/media/697054/tamu-cs-tr-2005-11-6.pdf
//---------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "opengl.hpp"
#include "geometry.hpp"

struct particle {
	cgra::vec3 pos;
	cgra::vec3 vel;
	cgra::vec3 acc;
	cgra::vec3 col;
	cgra::vec3 triangleIntersectionPos;
	int triangleIndex;
	bool inCollision;
};

class FuzzyObject {

	public:
		// Particle spawn position
		cgra::vec3 spawnPoint = cgra::vec3(0, 0, 0);

		// Constructors
		FuzzyObject(Geometry*);
		~FuzzyObject();

		void buildSystemIncrement();
		void buildSystem(bool);
		void addParticle();
		void renderSystem();

		int getParticleCount();
		void toggleParticleViewMode();

	private:
		// The 3D object the particle system represents
		Geometry* g_geometry;

		// Particle system fields
		std::vector<particle> particles;
		int particleLimit = 3000;
		std::vector<int> particlesForDeletion;

		// State fields
		bool buildFinished = false;

		// Stopping criteria
		bool firstPassFinished = false;
		float stabilityUpdates = 10;

		float manualShiftAmount = 0.001f;

		int collisionCount = 0;
		int lastCollisionCount = 0;

		int overfullStepCount = 0;
		int overfullThreshold = 10;
		bool overfull = false;

		int stabilityStepCount = 0;
		int stabilityThreshold = 30;
		bool wasStable = false;

		// Particle attributes
		GLuint p_displayList = 0;
		float p_velRange = 0.03f;
		float p_radius = 0.3f;
		float p_boundaryRadius = 0.05f;
		float p_mass = 100.0f;

		// LJ potential energy fields
		float e_strength = 0.005f;
		float e_lengthScale = 0.4f;
		float e_effectRange = pow(2.0f, 1.0f / 6.0f) * e_lengthScale;

		// Physics fields
		float meshCollisionFriction = 0.998f;
		float particleCollisionFriction = 0.998f;

		// Drawing properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;
		bool particleViewMode = true;

		cgra::vec3 maxFloatVector = cgra::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

		void setupDisplayList();
		bool stoppingCriteria();
		bool systemAtRest();
		void updateSystem();
		void applyParticleForces();
		void applyBoundaryForces();
		cgra::vec3 forceAtDistance(float, cgra::vec3);
		bool withinRange(cgra::vec3, cgra::vec3, float);
		void updateFacingTriangle(int);
};
