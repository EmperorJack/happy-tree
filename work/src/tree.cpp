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
	generateAttractionPointsVolumetric(500);
	root = generateTree();
}

branch* Tree::generateTree(){

	float d = param_branchLength;
	
	branch *root = new branch();
	branch *parent = root;
	branch *curNode = root;

	curNode->position = vec3(0,0,0);
	curNode->direction = vec3(0,1,0);
	curNode->length = (param_radiusOfInfluence + d < trunkHeight) ? trunkHeight : d;
	treeNodes.push_back(curNode);

	//Generate branches from attraction points
	// int prevSize = attractionPoints.size() + 1;
	while(attractionPoints.size() > 0){
		// cout << "treeSize " << treeNodes.size() << " attPoints " << attractionPoints.size() << endl;
		
		vector<vector<int>> closestSet = getAssociatedPoints();
		vector<branch *> toBeAdded;
		//Loop for all treeNodes
		for(int t=0; t<treeNodes.size(); t++){
			//Check if we want to branch
			if(closestSet[t].size() > 0){
				vec3 v = treeNodes[t]->position + (treeNodes[t]->direction * treeNodes[t]->length);
				vec3 newDir = vec3(0,0,0);

				for(int j=0; j<closestSet[t].size(); j++){
					int ind = closestSet[t][j];
					newDir += normalize(attractionPoints[ind] - v);
				}
				newDir = normalize(newDir);

				branch* newNode = new branch();
				newNode->position = v;
				newNode->direction = newDir;
				newNode->length = d;
				newNode->parent = treeNodes[t];
				treeNodes[t]->children.push_back(newNode);
				
				toBeAdded.push_back(newNode);
			}
		}
		treeNodes.insert(treeNodes.end(), toBeAdded.begin(), toBeAdded.end());
		cullAttractionPoints();
		// prevSize = attractionPoints.size();
	}

	setWidth(root);
	root->widthBase = root->widthTop;
	return root;
}

float Tree::setWidth(branch *b){
	float width = 0.0;

	//cout << "branch with children: " << b->children.size() << endl;

	for(int i=0; i<b->children.size(); i++){
		width += setWidth(b->children[i]);
	}

	width = (width == 0) ? 0.1 : width;

	b->widthTop = width;

	for(int i=0; i<b->children.size(); i++){
		(b->children[i])->widthBase = width;
	}

	return width;
}

vector<vector<int>> Tree::getAssociatedPoints(){
	vector<vector<int>> closestNodes;
	//init all the sets;
	for(int j=0; j<treeNodes.size(); j++){
		vector<int> sv = vector<int>();
		closestNodes.push_back(sv);
	}
	//Scan through all attraction points
	for(int i=0; i<attractionPoints.size(); i++){
		vec3 aPoint = attractionPoints[i];

		float minDist = distance(aPoint,treeNodes[0]->position + (treeNodes[0]->direction * treeNodes[0]->length));
		int closest = 0;

		for(int j=1; j<treeNodes.size(); j++){
			branch* t = treeNodes[j];
			vec3 p = t->position + (t->direction * t->length);
			float dist = distance(aPoint,p);
			if(dist <= minDist){
				closest = j;
				minDist = dist;
			}
		}

		//Only assign to the set if it is within the radius of influence
		if(minDist <= param_radiusOfInfluence){
			closestNodes[closest].push_back(i);
		}
	}
	return closestNodes;
}

void Tree::cullAttractionPoints(){
	int countRemoved = 0;

	for(int i=0; i<attractionPoints.size() - countRemoved;){
		vec3 aPoint = attractionPoints[i];

		bool toRemove = false;
		for(int j=0; j<treeNodes.size(); j++){
			branch* t = treeNodes[j];
			vec3 p = t->position + (t->direction * t->length);
			if(distance(aPoint,p) < param_killDistance){
				toRemove = true;
				break;
			}
		}
		if(toRemove){
			vec3 temp = attractionPoints[i];
			int ind = attractionPoints.size() - (1+countRemoved);
			attractionPoints[i] = attractionPoints[ind];
			attractionPoints[ind] = temp;

			countRemoved++;
		}else{
			i++;
		}
	}
	if(attractionPoints.size() <= countRemoved){
		attractionPoints.erase(attractionPoints.begin(), attractionPoints.end());
	}else{
		attractionPoints.erase(attractionPoints.end() - (1+countRemoved),attractionPoints.end());
	}
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
			float d = envelopeFunction(y-trunkHeight,theta);

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
	float uN = (2*u)/(treeHeight-trunkHeight);
	return 2*(pow(3,uN) - (uN*uN*uN));
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

	glPushMatrix();{
		vec3 basRot = b->basisRot;
		vec3 dir = b->direction;
		vec3 rot = b->rotation;
		/*
		//align local axis with world axis
		glRotatef(basRot.z,0,0,1);
		glRotatef(basRot.y,0,1,0);
		glRotatef(basRot.x,1,0,0);

		//rotate the body
		glRotatef(rot.z,0,0,1);
		glRotatef(rot.y,0,1,0);
		glRotatef(rot.x,1,0,0);

		//translate back to global axis
		glRotatef(-basRot.x,1,0,0);
		glRotatef(-basRot.y,0,1,0);
		glRotatef(-basRot.z,0,0,1);
		*/

		glPushMatrix();{
			float angle = acos(dot(normalize(dir),vec3(0,0,1)));
			vec3 axis = cross(normalize(dir),vec3(0,0,1));
			glRotatef(-degrees(angle),axis.x,axis.y,axis.z);

			cgraCylinder(b->widthBase, b->widthTop, b->length);

			//cgraSphere(b->widthBase);

		}glPopMatrix();
		
		
		vec3 offset = dir * b->length;
		glTranslatef(offset.x,offset.y,offset.z);

		for(branch* c : b->children){
			renderBranch(c);
		}
	}glPopMatrix();
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

void Tree::renderStick(){
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glTranslatef(m_position.x, m_position.y, m_position.z);

	//Actually draw the skeleton
	renderStick(root);

	// Clean up
	glPopMatrix();
}

void Tree::renderStick(branch *b){
	glPushMatrix();{
		glBegin(GL_LINE_STRIP);
		vec3 p1 = b->position;
		vec3 p2 = b->position + (b->direction * b->length);
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
		glEnd();

		for(branch* child : b->children){
			renderStick(child);
		}	
	}glPopMatrix();
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