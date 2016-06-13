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

struct fuzzyParticle {

	// Basic properties
	cgra::vec3 pos;
	cgra::vec3 vel;
	cgra::vec3 acc;
	cgra::vec3 col; // Colour

	// Collision properties
	cgra::vec3 triangleIntersectionPos;
	int triangleIndex;
	bool inCollision;

	// Neighbour relation properties
	int id;
	std::vector<int> neighbours;
};

class FuzzyObject {

	public:
		// Particle spawn position
		cgra::vec3 spawnPoint = cgra::vec3(0, 0, 0);

		// Constructors
		FuzzyObject(Geometry*);
		~FuzzyObject();

		// Methods for building the system
		void buildSystemIncrement();
		void buildSystem(bool);

		// Misc methods
		int getParticleCount();
		void setExampleSystemAttributes();

		// Methods for utilizing the built system
		bool finishedBuilding();
		std::vector<cgra::vec3> getSystem();
		void clearParticles();
		void scaleDensity(float);

		// Methods rendering the system
		void renderSystem();

	private:
		// The 3D object the particle system represents
		Geometry* g_geometry;

		// Particle system fields
		std::vector<fuzzyParticle> particles;
		int particleLimit = 3000;
		int minParticleCount = 10;
		std::vector<int> particlesForDeletion;
		int nextUniqueId = 0;

		// State fields
		bool buildFinished = false;

		// Stopping criteria
		bool firstPassFinished = false;
		float stabilityUpdates = 10;
		int collisionCount = 0;

		// Particle attributes
		GLuint p_displayList = 0;
		float p_velRange = 0.033f;
		float p_radius = 0.2f;
		float p_boundaryRadius = 0.16f;
		float p_mass = 100.0f;
		float p_spawnOffset = 0.018f;

		// LJ potential energy fields
		float e_strength = 0.01f;
		float e_lengthScale = 0.27f;
		float e_effectRange = pow(2.0f, 1.0f / 6.0f) * e_lengthScale;

		// Physics fields
		float meshCollisionFriction = 0.995f;
		float particleCollisionFriction = 0.995f;

		// Drawing properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 64.0f;

		cgra::vec3 maxFloatVector = cgra::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

		// Private methods for building the system
		void setupDisplayList();
		bool stoppingCriteria();
		bool systemAtRest();
		void addParticle();
		void updateBuildingSystem();
		void applyParticleForces();
		void applyBoundaryForces();
		cgra::vec3 forceAtDistance(float, cgra::vec3);
		bool withinRange(cgra::vec3, cgra::vec3, float);
		void updateFacingTriangle(int);
};
