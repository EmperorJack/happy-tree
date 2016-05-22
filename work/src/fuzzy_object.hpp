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

		// Constructors
		FuzzyObject(Geometry*);
		~FuzzyObject();

		void updateSystem();
		void renderSystem();

	private:
		// The 3D object the particle system represents
		Geometry* geometry;

		// Particle fields
		std::vector<particle> particles;
		int maxParticles = 2000;
		float pRadius = 0.2f;
		float pMass = 1.0f;

		// Material properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;

		void buildSystem();
};
