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
		void blowAway(cgra::vec3 direction);

		void resetParticles();

	private:
		
		// System fields
		std::vector<particle> particles;
		GLuint p_displayList = 0;

		// Particle fields
		float p_radius = 0.2f;
		float p_velRange = 0.03f;
		float p_maxVel = 0.15f;

		// Drawing properties
		cgra::vec4 diffuse = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		cgra::vec4 specular = cgra::vec4(0.8, 0.8, 0.8, 1.0);
		float shininess = 64.0f;

		// Animation fields
		int animationStep = 0;
		int animationLength = 100;
		cgra::vec4 currentColour = cgra::vec4(0.0, 0.0, 0.0, 1.0);
		cgra::vec4 startColour = cgra::vec4(1.0, 1.0, 1.0, 1.0);
		cgra::vec4 endColour = cgra::vec4(1.0, 0.2, 0.0, 0.0);

		void setupDisplayList();
};
