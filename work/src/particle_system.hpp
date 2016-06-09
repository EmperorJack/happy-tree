//---------------------------------------------------------------------------
// General Particle System
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

struct particle {
	cgra::vec3 original_pos;
	cgra::vec3 pos;
	cgra::vec3 vel;
	cgra::vec3 acc;
	cgra::vec3 col;
};

class ParticleSystem {

	public:

		// Constructors
		ParticleSystem(std::vector<cgra::vec3>);

		// Methods
		void update();
		void render();

		// Animation triggers
		void drop();
		void explode();

		void resetParticles();

	private:
		
		// System fields
		std::vector<particle> particles;
		GLuint p_displayList = 0;

		// Particle fields
		float p_radius = 0.2f;
		float p_velRange = 0.03f;

		// Drawing properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 128.0f;

		void setupDisplayList();
};
