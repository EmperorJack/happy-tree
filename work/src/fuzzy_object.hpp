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

		void buildSystem();
		void buildIncremental();
		void addParticle();
		void renderSystem();

		int getParticleCount();

	private:
		// The 3D object the particle system represents
		Geometry* geometry;
		std::vector<triangle> m_triangles;

		// Particle system fields
		std::vector<particle> particles;
		int particleLimit = 2000;

		// Particle attributes
		GLuint p_displayList = 0;
		float p_velRange = 0.01f;
		float p_radius = 0.1f;
		float p_boundaryRadius = 0.1f;
		float p_mass = 5.0f;

		// LJ potential energy fields
		float e_strength = 5.0f;
		float e_lengthScale = 0.1f;
		float e_effectRange = pow(2.0f, 1.0f / 6.0f) * e_lengthScale;

		float velocityReductionOnCollision = 0.9f;
		
		// Material properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;

		void setupDisplayList();
		void updateSystem();
		void applyParticleForces();
		void applyBoundaryForces();
		cgra::vec3 forceAtDistance(float dist, cgra::vec3);
		bool stoppingCriteria();
		bool systemAtRest();
};
