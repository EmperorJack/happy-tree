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
		void renderSystem();

	private:
		// The 3D object the particle system represents
		Geometry* geometry;

		// Particle system fields
		std::vector<particle> particles;
		int particleLimit = 5000;

		// Particle attributes
		float p_velRange = 0.01f;
		float p_radius = 0.1f;
		float p_mass = 1.0f;

		// LJ potential energy fields
		float strength = 1.0f;
		float lengthScale = 1.0f;
		float effectRange = pow(2.0f, 1.0f / 6.0f) * lengthScale;

		// Material properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;

		void addParticle();
		void updateSystem();
		bool stoppingCriteria();
		bool systemAtRest();
};
