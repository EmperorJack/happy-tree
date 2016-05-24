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

Tree::Tree(){
	//root = makeDummyTree(4);
	generateEnvelope(5.0, 20);
}

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

void Tree::generateEnvelope(float height, int steps){
	//vector<vec3> envelope;
	envelope.push_back(vec3(0,0,0));
	float step = height/steps;
	float y;
	for(int i = 1; i < steps; i++){
		y = i * step;
		for(float theta = 0; theta <= 360.0f; theta += 15.0f){
			float d = envelopeFunction(y,theta);
			envelope.push_back(vec3(d * sin(radians(theta)), y, d * cos(radians(theta))));
		}
	}
	envelope.push_back(vec3(0,height,0));

	
}

float Tree::envelopeFunction(float u, float theta){
	return ((u - 5) * u)/4;
}

void Tree::drawEnvelope(){
	glBegin(GL_LINE_STRIP);
	for(int i=0; i<envelope.size(); i++){
		vec3 p = envelope[i];
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}