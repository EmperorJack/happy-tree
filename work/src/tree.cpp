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
	setWindForce(vec3(20,0,20));
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
	root->baseWidth = root->topWidth;
	return root;
}

float Tree::setWidth(branch *b){
	float width = 0.0;

	//cout << "branch with children: " << b->children.size() << endl;

	for(int i=0; i<b->children.size(); i++){
		width += setWidth(b->children[i]);
	}

	width = (width == 0) ? 0.1 : width;

	b->topWidth = width;

	for(int i=0; i<b->children.size(); i++){
		(b->children[i])->baseWidth = width;
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

	//Actually draw the tree
	renderBranch(root);

	// Clean up
	glPopMatrix();
}

void Tree::renderBranch(branch *b) {
	if(b == NULL){
		return;
	}

	if(windEnabled){
		applyWind(b);
	}
	glPushMatrix();

		//only draw branch info if it has a length
		if(b->length > 0){
		
			vec3 rot = b->basisRot;
			glRotatef(rot.z, 0, 0, 1);
			glRotatef(rot.y, 0, 1, 0);
			glRotatef(rot.x, 1, 0, 0);

			cout << b->name << endl;
			cout << "Branch Rotation X: " <<  b->rotation.x << endl;
			//cout << "Branch Rotation Y: " <<  b->rotation.y << endl;
			cout << "Branch Rotation Z: " <<  b->rotation.z << endl;
			cout << endl;

			glRotatef(b->rotation.z, 0, 0, 1);
			//glRotatef(b->rotation.y, 0, 1, 0);
			glRotatef(b->rotation.x, 1, 0, 0);
	
			glRotatef(-rot.x, 1, 0, 0);
			glRotatef(-rot.y, 0, 1, 0);
			glRotatef(-rot.z, 0, 0, 1);

		glPushMatrix();{
			//draw the joint of this branch
			//drawJoint(b);
			drawBranch(b);
		}glPopMatrix();
		
		
		vec3 offset = b->direction * b->length;
		glTranslatef(offset.x,offset.y,offset.z);

		for(branch* c : b->children){
			renderBranch(c);
		}
	}glPopMatrix();
}

void Tree::drawJoint(branch* b){
	glPushMatrix();
		//colour cyan
		glColor3f(0,1,1);
		cgraSphere(b->baseWidth);
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
		cgraCylinder(b->baseWidth, b->topWidth, b->length);
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
	b->baseWidth = width * numBranches;
	b->topWidth = (width * (numBranches - 1)) + 0.01f;
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
			c->baseWidth = width * (numBranches-1);
			c->topWidth = 0.01f;
			c->basisRot = vec3(0,0,0);


			b->children.push_back(c);
		}

		b->children.push_back(makeDummyTree(numBranches - 1));

		
	}
	return b;
}

void Tree::setWindForce(vec3 wind){
	windForce = wind;
}

float Tree::calculatePressure(branch* branch, float force){
	float t = branch->baseWidth - branch-> topWidth;

	float a = 1.0f; //change to a small number derived from the current angle of the branch
	//float b = math::random(0.0f,0.1f); //change to random value to make different branches be at different point of sine equation

	float oscillation = (time+branch->offset);

	//remove *degrees if calculation should be in radians
	float degrees =  (float)(math::pi());
	degrees = degrees / 180.0f;

	float pressure = force * (1 + a * sin(oscillation) );

	return pressure;
}

float Tree::springConstant(branch* branch){
	float k = (elasticity*branch->baseWidth*pow(branch->baseWidth-branch->topWidth, 3));
	k = k/ (4*pow(branch->length, 3));

	return k;
}

float Tree::displacement(branch* branch, float pressure){
	float spring = springConstant(branch);

	return pressure/spring;
}

void Tree::applyWind(branch* b){
	time += 0.000008f;

	float displacementX = displacement(b, calculatePressure(b, (windForce.x)));
	//float displacementY = displacement(b, calculatePressure(b, (windForce.y)));
	float displacementZ = displacement(b, calculatePressure(b, (windForce.z)));

	cout << "length " << b->length << endl;
	cout << "Displacement - x: " << displacementX /*<< "  y: " << displacementY */<< "  z: " << displacementZ << endl;
	
	float degrees =  ( (float)(math::pi()) ) / 180.0f;

	int len = b->length;

	if(len == 0){
		len = 0.00001f;
	}

	float motionAngleX = asin(displacementX/float(len)); //* degrees;
	//float motionAngleY = asin(displacementY/float(len)); //* degrees;
	float motionAngleZ = asin(displacementZ/float(len));// * degrees;

	cout << "Motion Angle - x: " << motionAngleX /*<< "  y: " << motionAngleY */<< "  z: " << motionAngleZ << endl;
	cout << asin(5) << endl;
	cout << endl;

	b->rotation.x = motionAngleX;
	//b->rotation.y = motionAngleY;
	b->rotation.z = motionAngleZ;

	b->rotation.x = displacementX;
	//b->rotation.y = displacementY;
	b->rotation.z = displacementZ;
}

void Tree::toggleWind(){
	windEnabled = ! windEnabled;
}