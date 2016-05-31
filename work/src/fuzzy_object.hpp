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
		Geometry* geometry;
		std::vector<triangle> m_triangles;

		// Particle system fields
		std::vector<particle> particles;
		int particleLimit = 5000;

		// Stopping criteria
		bool finished = false;
		float manualShiftAmount = 0.001f;

		int collisionCount = 0;
		int lastCollisionCount = 0;

		int overfullStepCount = 0;
		int overfullThreshold = 10;
		bool overfull = false;

		int stabilityStepCount = 0;
		int stabilityThreshold = 80;
		bool wasStable = false;

		// Particle attributes
		GLuint p_displayList = 0;
		float p_velRange = 0.02f;
		float p_radius = 0.05f;
		float p_boundaryRadius = 0.05f;
		float p_mass = 100.0f;

		// LJ potential energy fields
		float e_strength = 0.001f;
		float e_lengthScale = 0.3f;
		float e_effectRange = pow(2.0f, 1.0f / 6.0f) * e_lengthScale;

		// Physics fields
		float meshCollisionFriction = 0.6f;
		float particleCollisionFriction = 0.95f;

		// Drawing properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;
		bool particleViewMode = true;

		void setupDisplayList();
		bool stoppingCriteria();
		bool systemAtRest();
		void updateSystem();
		void applyParticleForces();
		void applyBoundaryForces();
		cgra::vec3 forceAtDistance(float dist, cgra::vec3);
};
