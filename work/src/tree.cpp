#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "cgra_geometry.hpp"
#include "cgra_math.hpp"
#include "opengl.hpp"
#include "tree.hpp"


using namespace std;
using namespace cgra;

float treeHeight;
float yStep;

Tree::Tree(){
	treeHeight = 8.0f;
	generateEnvelope(20);
}

void Tree::generateEnvelope(int steps){
	vector<vector<vec3>> env;
	// vector base;
	// base.push_back(vec3(0,0,0));
	// env.push_back(base);

	yStep = treeHeight/steps;
	float y;

	for(int i = 0; i <= steps; i++){
		vector<vec3> layer;
		y = i * yStep;
		for(float theta = 0; theta <= 360.0f; theta += 15.0f){
			float d = envelopeFunction(y,theta,treeHeight);
			layer.push_back(vec3(d * sin(radians(theta)), y, d * cos(radians(theta))));
		}
		env.push_back(layer);
	}

	envelope = env;
}

float Tree::envelopeFunction(float u, float theta, float range){
	return u < (range * 1.0f/3.0f) ? 1 : 3;
}

void Tree::drawEnvelope(){
	for(int i=0; i<envelope.size(); i++){
		vector<vec3> layer = envelope[i];

		glBegin(GL_LINE_STRIP);
		for(int j=0; j<layer.size(); j++){
			vec3 p = layer[j];
			glVertex3f(p.x, p.y, p.z);
		}
		glEnd();
	}
}

// UNNECESSARY METHOD
// Only used so that there is a model to work with
branch* Tree::makeDummyTree(int numBranches){
	branch b;
	b.direction = vec3(0,1,0);
	b.length = 2.0f;
	b.widthBase = 0.5f;
	b.widthTop = 0.4f;
	if(numBranches > 1){
		b.children.push_back(makeDummyTree(numBranches - 1));
	}
	return &b;
}