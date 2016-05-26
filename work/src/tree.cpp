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
	treeHeight = 15.0f;
	trunkHeight = 5.0f;
	generateEnvelope(20);
	generateAttractionPointsVolumetric(1000);
	//cout << inEnvelope(vec3(2.4,7.5f,2.4)) << endl;
}

void Tree::generateAttractionPoints(int numPoints){
	if(numPoints == 0) return;

	vector<vec3> points;
	while(points.size() < numPoints){
		//Calculate random height/rotation
		float y = math::random(trunkHeight,treeHeight);
		float theta = math::random(0.0f,360.0f);

		// Calculate max distance
		float d = envelopeFunction(y,theta);

		// Calculate distance away from central axis
		float r = math::random(0.0f,d);

		// Convert from rotation/distance to x,z
		points.push_back(vec3(r * sin(radians(theta)), y, r * cos(radians(theta))));
	}
	attractionPoints = points;
}

void Tree::generateAttractionPointsVolumetric(int numPoints){
	if(numPoints == 0) return;

	vector<vec3> points;
	while(points.size() < numPoints){
		float x = math::random(minX,maxX);
		float y = math::random(trunkHeight,treeHeight);
		float z = math::random(minZ,maxZ);

		vec3 point = vec3(x,y,z);

		if(inEnvelope(point)){
			points.push_back(point);
		}
	}
	attractionPoints = points;
}






//------------------------------------------------//
//   Envelope Functions                           //
//------------------------------------------------//
void Tree::generateEnvelope(int steps){
	vector<vector<vec3>> env;

	yStep = (treeHeight - trunkHeight)/steps;
	float y;

	for(int i = 0; i <= steps; i++){
		vector<vec3> layer;
		y = (i * yStep) + trunkHeight;
		for(float theta = 0; theta <= 360.0f; theta += 15.0f){
			float d = envelopeFunction(y,theta);

			float x = d * sin(radians(theta));
			float z = d * cos(radians(theta));
			
			//Assign bounding values for volumetric filling
			minZ = z < minZ ? z : minZ;
			maxZ = z > maxZ ? z : maxZ;
			minX = z < minX ? z : minX;
			maxX = z > maxX ? z : maxX;

			layer.push_back(vec3(x, y, z));
		}
		env.push_back(layer);
	}
	envelope = env;
}

bool Tree::inEnvelope(vec3 point){
	float x = point.x;
	float y = point.y;
	float z = point.z;
	//Make sure is within y bounds of envelope
	if(y < trunkHeight || y > treeHeight){
		return false;
	}
	//Calculate the level of the point;
	int yInd1 = int((y - trunkHeight)/yStep);
	int yInd2 = yInd1 + 1;
	// Ratio between y of points at yInd1/2
	float deltaY = (y - ((yInd1 * yStep) + trunkHeight))/yStep;

	vector<vec3> layer1 = envelope[yInd1];
	vector<vec3> layer2 = envelope[yInd2];

	float radius = distance(vec3(0,y,0), point);
	float theta = degrees(asin(x/radius));
	theta = theta < 0.0f ? 360.0f + theta : theta;

	// Calculate indicies based off rotation
	int xzInd1 = int(theta/thetaStep);
	int xzInd2 = xzInd1 + 1;
	// Ratio between the rotations of points at xzInd1/2
	float deltaT = (theta - (xzInd1 * thetaStep))/thetaStep;

	// Calculate the points between the two points with the y ratio
	vec3 xzP1 = layer1[xzInd1] + (deltaY * (layer2[xzInd1] - layer1[xzInd1]));
	vec3 xzP2 = layer1[xzInd2] + (deltaY * (layer2[xzInd2] - layer1[xzInd2]));

	// Calculate the point between the two points with the rotation ratio
	vec3 maxP = xzP1 + (deltaT * (xzP2 - xzP1));
	float maxRadius = distance(vec3(0,y,0), maxP);

	return radius <= maxRadius;
}

float Tree::envelopeFunction(float u, float theta){
	return 3;
}






//------------------------------------------------//
//   Rendering Functions                          //
//------------------------------------------------//

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

void Tree::renderAttractionPoints(){
	for(int i=0; i< attractionPoints.size(); i++){
		glPushMatrix();
		vec3 p = attractionPoints[i];
		glTranslatef(p.x,p.y,p.z);
		cgraSphere(0.1);
		glPopMatrix();
	}
}

void Tree::renderTree() {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glTranslatef(m_position.x, m_position.y, m_position.z);

	//Actually draw the skeleton
	renderBranch(root);

	// Clean up
	glPopMatrix();
}

void Tree::renderBranch(branch *b) {
	if(b == NULL){
		return;
	}

	glPushMatrix();

		//only draw branch info if it has a length
		if(b->length > 0){
		
			vec3 rot = b->basisRot;
			glRotatef(rot.z, 0, 0, 1);
			glRotatef(rot.y, 0, 1, 0);
			glRotatef(rot.x, 1, 0, 0);
		
			glRotatef(b->rotation.x, 1, 0, 0);
			glRotatef(b->rotation.y, 0, 1, 0);
			glRotatef(b->rotation.z, 0, 0, 1);
	
			glRotatef(-rot.x, 1, 0, 0);
			glRotatef(-rot.y, 0, 1, 0);
			glRotatef(-rot.z, 0, 0, 1);

			//glDisable(GL_LIGHTING);
			//draw the axes of this branch
			//drawAxis(b);

			//draw the joint of this branch
			//drawJoint();

			//glEnable(GL_LIGHTING);

			//draw the branch itself
			drawBranch(b);

			vec3 dir = b->direction;
			//translate to the end of the branch based off length and direction
			glTranslatef(dir.x*b->length, dir.y*b->length, dir.z*b->length);

		}
		//loop through all child branches and render them too
		for(branch* child : b->children){
			renderBranch(child);
		}	
	glPopMatrix();
}

void Tree::drawJoint(){
	glPushMatrix();
		//colour cyan
		glColor3f(0,1,1);
		cgraSphere(1.2*width);
	glPopMatrix();
}

void Tree::drawBranch(branch* b){
	vec3 dir = b->direction;
	vec3 vec = vec3(0,0,1);
	vec3 norm = normalize(dir);
	float dotProd = dot(norm, vec);
	float angle = acos(dotProd);
	vec3 crossProd = cross(dir, vec);
	
	glPushMatrix();
	//colour grey
		glColor3f(1,1,1);
		glRotatef(-degrees(angle), crossProd.x, crossProd.y, crossProd.z);
		cgraCylinder(b->widthBase, b->widthTop, b->length);
	glPopMatrix();
}

void Tree::drawAxis(branch* b){
	//X-axes
	glPushMatrix();
		glColor3f(1,0,0);
		glRotatef(90,0, 1,0);
		cgraCylinder(0.3*width, 0.3*width, 4.0*width);
		glTranslatef(0, 0, 4.0*width);
		cgraCone(1.2*width, 1.2*width);
	glPopMatrix();

	//Y-axes
	glPushMatrix();
		glColor3f(0,1,0);
		glRotatef(-90,1,0,0);
		cgraCylinder(0.3*width, 0.3*width, 4.0*width);
		glTranslatef(0,0, 4.0*width);
		cgraCone(1.2*width, 1.2*width);
	glPopMatrix();

	//Z-axes
	glPushMatrix();
		glColor3f(0,0,1);
		cgraCylinder(0.3*width, 0.3*width, 4.0*width);		
		glTranslatef(0,0,4.0*width);
		cgraCone(1.2*width, 1.2*width);	
	glPopMatrix();
}






//------------------------------------------------//
//   Miscellaneous Functions                      //
//------------------------------------------------//
void Tree::setPosition(vec3 position) {
	m_position = position;
}

// UNNECESSARY METHOD
// Only used so that there is a model to work with
branch* Tree::makeDummyTree(int numBranches){
	branch* b = new branch();
	b->direction = vec3(0,1,0);
	b->length = length;
	b->widthBase = width * numBranches;
	b->widthTop = (width * (numBranches - 1)) + 0.01f;
	b->basisRot = vec3(0,0,0);
	if(numBranches > 1){

		for (int i = 0; i < 4; i++){
			branch* c = new branch();

			if(i == 0){
				c->direction = vec3(1,0,0);
			}else if(i == 1){
				c->direction = vec3(-1,0,0);
			}else if(i == 2){
				c->direction = vec3(0,0,1);
			}else if(i == 3){
				c->direction = vec3(0,0,-1);
			}

			c->length = length;
			c->widthBase = width * (numBranches-1);
			c->widthTop = 0.01f;
			c->basisRot = vec3(0,0,0);


			b->children.push_back(c);
		}

		b->children.push_back(makeDummyTree(numBranches - 1));

		
	}
	return b;
}