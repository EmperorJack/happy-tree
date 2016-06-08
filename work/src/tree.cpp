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
	treeHeight = 20.0f;
	trunkHeight = 2.0f;

	param_branchLength = 2.0f;
	param_radiusOfInfluence = 8 * param_branchLength;
	param_killDistance = param_branchLength;
	param_branchTipWidth = 0.06;
	param_branchMinWidth = 0.08;

	generateEnvelope(20);
	generateAttractionPointsVolumetric(200);

	generatedTreeRoot = generateTree();
	generateGeometry(generatedTreeRoot);
	dummyTreeRoot = makeDummyTree(4); // make dummy tree to work with

	if(dummyTree){
		root = dummyTreeRoot;
	} else {
		root = generatedTreeRoot;
	}
}

branch* Tree::generateTree(){
	float d = param_branchLength;

	branch *root = new branch();
	branch *parent = root;
	branch *curNode = root;
	curNode->position = vec3(0,0,0);
	curNode->direction = vec3(0,1,0);
	curNode->length = trunkHeight;
	treeNodes.push_back(curNode);

	// while(curNode->position.y + d < trunkHeight){
	// 	curNode = new branch();
	// 	curNode->position = parent->position + (parent->direction * parent->length);
	// 	curNode->direction = vec3(0,1,0);
	// 	curNode->length = d;
	// 	curNode->parent = parent;
	// 	parent->children.push_back(curNode);
	// 	treeNodes.push_back(curNode);

	// 	parent = curNode;
	// }

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
				newDir = normalize(newDir + vec3(0,-0.1,0));

				branch* newNode = new branch();
				newNode->position = v;
				newNode->direction = newDir;
				newNode->length = d;
				newNode->parent = treeNodes[t];
				newNode->offset = math::random(0.0f,1.0f);

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
	float maxW = param_branchTipWidth;

	//cout << "branch with children: " << b->children.size() << endl;

	for(int i=0; i<b->children.size(); i++){
		float cw = setWidth(b->children[i]);
		width += pow(cw, 2);
		maxW = (cw > maxW) ? cw : maxW;
	}

	width = (width == 0) ? param_branchMinWidth : sqrt(width);

	b->topWidth = maxW;
	b->baseWidth = width;



	// for(int i=0; i<b->children.size(); i++){
	// 	(b->children[i])->baseWidth = width;
	// }

	return width;
}

void Tree::generateGeometry(branch *b) {
	b->jointModel = generateSphereGeometry(b->baseWidth);
	b->branchModel = generateCylinderGeometry(b->baseWidth, b->topWidth, b->length);

	b->jointModel->setMaterial(vec4(0.2, 0.2, 0.2, 1.0), vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, vec4(0.0, 0.0, 0.0, 1.0));
	b->branchModel->setMaterial(vec4(0.2, 0.2, 0.2, 1.0), vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, vec4(0.0, 0.0, 0.0, 1.0));

	for (branch* c : b->children) {
		generateGeometry(c);
	}
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
		for(float theta = 0; theta <= 360.0f; theta += thetaStep){
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
	float uN = u/(treeHeight-trunkHeight);
	// return 6*(pow(3,2*uN) - (8*uN*uN*uN));
	// return (1.0f - uN) * 8;
	return -100 * (uN * uN * (uN - 1));
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

/* public method for drawing the tree to the screen.
	draws the tree by calling renderbranch() on the root node.
*/
void Tree::renderTree() {
	//glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//makes sure the tree is drawn at its set position
	glTranslatef(m_position.x, m_position.y, m_position.z);

	//Actually draw the tree
	renderBranch(root);

	//increment wind "time"
	time += timeIncrement;

	// Clean up
	glPopMatrix();
}

/* performs the logic for drawing any given branch at its position and rotation.
	then recursively calls renderBranch() on all of its child branches.
*/
void Tree::renderBranch(branch *b, int depth) {
	if(b == NULL){
		return;
	}
	//togglable for starting and stopping the wind being applied
	if(windEnabled){
		applyWind(b);
	}

	static int i_k =0;
	glPushMatrix();
		//only draw branch info if it has a length
		if(b->length > 0){

			vec3 rot = b->basisRot;
			glRotatef(rot.z, 0, 0, 1);
			glRotatef(rot.y, 0, 1, 0);
			glRotatef(rot.x, 1, 0, 0);

			//debug info
			// cout << b->name << endl;
			// cout << "Branch Rotation X: " <<  b->rotation.x << endl;
			// cout << "Branch Rotation Z: " <<  b->rotation.z << endl;
			// cout << endl;

			//perform rotation as updated by wind
			glRotatef(b->rotation.z, 0, 0, 1);
			glRotatef(b->rotation.x, 1, 0, 0);

			glRotatef(-rot.x, 1, 0, 0);
			glRotatef(-rot.y, 0, 1, 0);
			glRotatef(-rot.z, 0, 0, 1);

			//draw the joint of this branch
			drawJoint(b);

			drawBranch(b);

			//translate to the end of the branch based off length and direction
			vec3 offset = b->direction * b->length;
			glTranslatef(offset.x,offset.y,offset.z);
		}

		// if (depth < 20 && i_k < 1000) {
		// 	i_k++;
		// 	for (int i = 0; i < depth; ++i) cout << " ";
		// 	cout << "branch :: pos " << b->position << " : post " << pos << " : offset " << offset << endl;
		// }
		//loop through all child branches and render them too

		for(branch* c : b->children){
			renderBranch(c, depth+1);
		}

	glPopMatrix();
}

void Tree::renderStick(){
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	// glTranslatef(m_position.x, m_position.y, m_position.z);

	//Actually draw the skeleton
	renderStick(root);

	// Clean up
	glPopMatrix();
}

void Tree::renderStick(branch *b, int depth){
	glPushMatrix();{
		int n = depth * 15;
		int cR = (n > 255) ? 255 : n;
		int cB = (n >= 255) ? ( n > 510 ? 255 : (n - 255) ) : 0;
		int cG = (n >= 510) ? n - 510 : 0;

		glColor3f(cR / 255.0f,cB / 255.0f,cG / 255.0f);
		glBegin(GL_LINES);
		vec3 p1 = b->position;
		vec3 p2 = b->position + (b->direction * b->length);
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
		glEnd();

		// vec3 dir = b->direction;
		// vec3 offset = dir * b->length;
		// glBegin(GL_LINES);
		// 	glVertex3f(0.0f,0.0f,0.0f);
		// 	glVertex3f(offset.x,offset.y,offset.z);
		// glEnd();
		// glTranslatef(offset.x,offset.y,offset.z);

		for(branch* child : b->children){
			renderStick(child, depth+1);
		}
	}glPopMatrix();
}

/* draws a joint at the base of every branch the size of the width at the base of the branch
	this prevents a tree breaking visual issue when rotating branches.
*/
void Tree::drawJoint(branch* b){
	glPushMatrix();
		//cgraSphere(b->baseWidth);
		b->jointModel->renderGeometry();
	glPopMatrix();
}

/* draws the branch to the screen
*/
void Tree::drawBranch(branch* b){
	vec3 norm = normalize(b->direction);
	float dotProd = dot(norm, vec3(0,0,1));

	float angle = acos(dotProd); // the angle to rotate by
	vec3 crossProd = cross(b->direction, vec3(0,0,1));

	glPushMatrix();
		glRotatef(-degrees(angle), crossProd.x, crossProd.y, crossProd.z);
		//cgraCylinder(b->baseWidth, b->topWidth, b->length);
		b->branchModel->renderGeometry();
	glPopMatrix();
}

/*
	Calculates the pressure the wind will apply to a given branch
	force is the float value of the wind in the windforce vector for a given axis (x or z)
	dir is an int to let us know which axis to calculate the force for (0 == x, 2 == z)
*/
float Tree::calculatePressure(branch* branch, float force, int dir){
	float a = windCoefficent; //change to a small number derived from the current angle of the branch

	//attempt at making the small value use the current angle of the branch
	// if (dir == 'x'){ //x axis
	// 	a = sin(branch->rotation.z);
	// } else if (dir == 'z'){ //z axis
	// 	a = sin(branch->rotation.x);
	// }

	//oscillation is plugged into a sine function.
	//time is increased steadily to make the effect follow an oscilation pattern - global scope
	//branch offset is a random value assigned to each branch so they are at a different point in the oscillation
	float oscillation = (time + branch->offset);
	//oscillation = (oscillation - floor(oscillation)) - 0.5f;

	//mulitply a radians value by degrees variable to convert it from radians to degrees
	//not sure if this is needed...
	float degrees = 180.0f / ((float)math::pi()) ;

	//pressure is the final return value
	float pressure = force * (1 + a * sin(oscillation) );
	//float pressure = sin(oscillation);

	return pressure;
}


/*
	A spring value for a branch based on its thickness and length.
	Taken from a reserch paper
*/
float Tree::springConstant(branch* branch){
	float thickness = (branch->baseWidth+branch->topWidth)/2.0f;

	float k = (elasticity * branch->baseWidth *	pow(thickness, 3));

	// cout << "Spring top: " << k << endl;
	// cout << "Spring bot: " << (4 * pow( branch->length, 3)) << endl;

	k = k / (4 * pow( branch->length, 3));

	return k;
}

/*
	the central method for applying wind force to a branch.
	calculates the displacement value for the branch based on the wind then
	stores the value to rotate it by
*/
void Tree::applyWind(branch* b){
	//increment time (this value is currently has no meaning, just seems to fit at an ok speed)

	//calculates the pressure value for each axis
	float pressureX = calculatePressure(b, desiredWindForce.x, 'x');
	float pressureZ = calculatePressure(b, desiredWindForce.z, 'z');

	//debug info
	// cout << "Name: " << b->name << endl;
	// cout << "Pressure X: " << pressureX << endl;
	// cout << "Pressure Z: " << pressureZ << endl;

	//the spring value of this branch
	float spring = springConstant(b);
	// cout << "Spring Value: " << spring << endl;

	//make sure no division of 0 is occuring
	int len = b->length;
	if(len == 0){
		len = 0.00001f;
	}
	if(spring == 0){
		spring = 0.00001f;
	}

	//calculates the displacement value for each axis
	float displacementX = pressureX / spring / float(len);
	float displacementZ = pressureZ / spring / float(len);

	//debug info
	// cout << "length " << b->length << endl;
	// cout << "Displacement - x: " << displacementX << "  z: " << displacementZ << endl;

	//clamp values to be within [-1,1]
	if (displacementX > 1){
		displacementX = 1.0f;
	} else if (displacementX < -1){
		displacementX = -1.0f;
	}

	if (displacementZ > 1){
		displacementZ = 1.0f;
	} else if (displacementZ < -1){
		displacementZ = -1.0f;
	}

	float motionAngleX = asin(displacementX);
	float motionAngleZ = asin(displacementZ);

	if(motionAngleX > b->maxX){
		b->maxX = motionAngleX;
	} else if (motionAngleX < b->minX){
		b->minX = motionAngleX;
	}
	if(motionAngleZ > b->maxZ){
		b->maxZ = motionAngleZ;
	} else if (motionAngleZ < b->minZ){
		b->minZ = motionAngleZ;
	}

	// cout << "Min Angle - x: " << b->minX << "  z: " << b->minX << endl;
	// cout << "Max Angle - x: " << b->maxZ << "  z: " << b->maxZ << endl;

	// cout << "Motion Angle - x: " << motionAngleX << "  z: " << motionAngleZ << endl;

	//mulitply a radians value by degrees variable to convert it from radians to degrees
	//not sure if this is needed...
	float degrees = 180.0f / ((float)math::pi()) ;

	b->rotation.z = motionAngleX * degrees;
	b->rotation.x = motionAngleZ * degrees;

	if (b->rotation.z > 20){
		b->rotation.z = 20.0f;
	} else if (b->rotation.z < -20){
		b->rotation.z = -20.0f;
	}

	if (b->rotation.x > 20){
		b->rotation.x = 20.0f;
	} else if (b->rotation.x < -20){
		b->rotation.x = -20.0f;
	}
	//temporarily just rotating by displacement value because motionAngle is NaN
	//attempt to restrict the rotation by converting it to degrees, and then limit it to 20degrees
	//b->rotation.x = ((displacementX*degrees) / 180 ) * 20;
	//b->rotation.z = ((displacementZ*degrees) / 180 ) * 20;

	//debug info
	// cout << "Final Angle - x: " << b->rotation.x << "  z: " << b->rotation.z << endl;
	// cout << endl;

	// if(b->rotation.x > 20){
	// 	b->rotation.x = 20;
	// } else if(b->rotation.x < -20){
	// 	b->rotation.x = -20;
	// }

	// if(b->rotation.z > 20){
	// 	b->rotation.z = 20;
	// } else if(b->rotation.z < -20){
	// 	b->rotation.z = -20;
	// }

}

//------------------------------------------------//
//   Miscellaneous Functions                      //
//------------------------------------------------//
/*
	sets the position we want to draw the tree at
*/
void Tree::setPosition(vec3 position) {
	m_position = position;
}

/*
	public method for toggling wind off and on
*/
void Tree::toggleWind(){
	windEnabled = ! windEnabled;
}

/*
	public method for toggling between dummy tree model and randomly generated tree
*/
void Tree::toggleTreeType(){
	dummyTree = ! dummyTree;
	if(dummyTree){
		root = dummyTreeRoot;
	} else {
		root = generatedTreeRoot;
	}
}

/*
	public method for generating a new tree
*/
void Tree::generateNewTree(){
	treeHeight = 20.0f;
	trunkHeight = 2.0f;

	param_branchLength = 2.0f;
	param_radiusOfInfluence = 8 * param_branchLength;
	param_killDistance = param_branchLength;
	param_branchTipWidth = 0.06;
	param_branchMinWidth = 0.08;

	generateEnvelope(20);
	generateAttractionPointsVolumetric(200);

	generatedTreeRoot = generateTree();

	if(dummyTree){
		root = dummyTreeRoot;
	} else {
		root = generatedTreeRoot;
	}
}

/*
	sets the wind force
*/
void Tree::setWindForce(vec3 wind){
	desiredWindForce = wind;
}

void Tree::adjustWind(int axis, int dir){
	float wIncrease = 0.00005f;
	float aIncrease = 0.1f;
	float tIncrease = 0.002f;

	if (axis == 'x'){
		if (dir == 1){
			desiredWindForce.x += wIncrease;
		} else if (dir == -1){
			desiredWindForce.x -= wIncrease;
		}
	} else if (axis == 'z'){
		if (dir == 1){
			desiredWindForce.z += wIncrease;
		} else if (dir == -1){
			desiredWindForce.z -= wIncrease;
		}
	}  else if (axis == 'a'){
		if (dir == 1){
			windCoefficent += aIncrease;
		} else if (dir == -1){
			windCoefficent -= aIncrease;
		}
	}  else if (axis == 't'){
		if (dir == 1){
			timeIncrement += tIncrease;
		} else if (dir == -1){
			timeIncrement -= tIncrease;
		}
	}
}

vector<Geometry*> Tree::getGeometries() {
	vector<Geometry*> geometries;

	getBranchGeometry(root, &geometries);

	return geometries;
}

void Tree::getBranchGeometry(branch* b, vector<Geometry*>* geometries) {
	geometries->push_back(b->jointModel);
	geometries->push_back(b->branchModel);

	for (branch* c : b->children) {
		getBranchGeometry(b, geometries);
	}
}

/* Builds a test tree to work with for simulating wind animation.
	Trees is 'numBranches' segments tall, with 4 branches inbetween each segment.
*/
branch* Tree::makeDummyTree(int numBranches){
	//hardcoded values for this dummy tree
	float width = 0.1f;
	float length = 5.0f;

	branch* b = new branch();
	b->name = "trunk"+to_string(numBranches);
	b->direction = vec3(0,1,0);
	b->offset = math::random(0.0f,1.0f);
	b->length = length;
	b->baseWidth = width * numBranches;
	b->topWidth = (width * (numBranches - 1));
	if(numBranches == 1){
		b->topWidth = 0.0001f;//(width /2);
	}

	b->basisRot = vec3(0,0,0);
	if(numBranches > 1){

		for (int i = 0; i < 4; i++){
			branch* c = new branch();

			c->name = "branch"+to_string(i+1)+" trunk"+to_string(numBranches);

			if(i == 0){
				c->direction = vec3(1,0.3,0);
			}else if(i == 1){
				c->direction = vec3(-1,0.3,0);
			}else if(i == 2){
				c->direction = vec3(0,0.3,1);
			}else if(i == 3){
				c->direction = vec3(0,0.3,-1);
			}

			c->length = length/2 * (numBranches-1);
			c->baseWidth = width * (numBranches-1);
			c->topWidth = width/2;
			c->basisRot = vec3(0,0,0);


			b->children.push_back(c);
		}

		b->children.push_back(makeDummyTree(numBranches - 1));


	}
	return b;
}
