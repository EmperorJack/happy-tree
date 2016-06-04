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
	root = makeDummyTree(4); // make dummy tree to work with
	setWindForce(vec3(20,0,20)); //set wind 
}

/* public method for drawing the tree to the screen.
	draws the tree by calling renderbranch() on the root node.
*/
void Tree::renderTree() {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//makes sure the tree is drawn at its set position
	glTranslatef(m_position.x, m_position.y, m_position.z);

	//Actually draw the tree
	renderBranch(root);

	// Clean up
	glPopMatrix();
}

/* performs the logic for drawing any given branch at its position and rotation.
	then recursively calls renderBranch() on all of its child branches.
*/
void Tree::renderBranch(branch *b) {
	if(b == NULL){
		return;
	}
	//togglable for starting and stopping the wind being applied
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

			//debug info
			cout << b->name << endl;
			cout << "Branch Rotation X: " <<  b->rotation.x << endl;
			cout << "Branch Rotation Z: " <<  b->rotation.z << endl;
			cout << endl;

			//perform rotation as updated by wind
			glRotatef(b->rotation.z, 0, 0, 1);
			glRotatef(b->rotation.x, 1, 0, 0);
	
			glRotatef(-rot.x, 1, 0, 0);
			glRotatef(-rot.y, 0, 1, 0);
			glRotatef(-rot.z, 0, 0, 1);

			//draw the joint of this branch
			drawJoint(b);

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

/* draws a joint at the base of every branch the size of the width at the base of the branch
	this prevents a tree breaking visual issue when rotating branches.
*/
void Tree::drawJoint(branch* b){
	glPushMatrix();
		cgraSphere(b->baseWidth);
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
		cgraCylinder(b->baseWidth, b->topWidth, b->length);
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
	if (dir == 'x'){ //x axis
		//a = branch->rotation.x;
	} else if (dir == 'z'){ //z axis
		//a = branch->rotation.z;
	}

	//oscillation is plugged into a sine function. 
	//time is increased steadily to make the effect follow an oscilation pattern - global scope
	//branch offset is a random value assigned to each branch so they are at a different point in the oscillation
	float oscillation = (time + branch->offset);

	//mulitply a radians value by degrees variable to convert it from radians to degrees
	//not sure if this is needed...
	float degrees =  ((float)(math::pi()))/180.0f;

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
	float thickness = branch->baseWidth-branch->topWidth;

	float k = (elasticity * branch->baseWidth *	pow(thickness, 3));
	k = k/ (4 * pow( branch->length, 3));

	return k;
}

/*
	the central method for applying wind force to a branch.
	calculates the displacement value for the branch based on the wind then
	stores the value to rotate it by
*/
void Tree::applyWind(branch* b){
	//increment time (this value is currently has no meaning, just seems to fit at an ok speed)
	time += 0.000008f;

	//calculates the pressure value for each axis	
	float pressureX = calculatePressure(b, windForce.x, 'x');
	float pressureZ = calculatePressure(b, windForce.z, 'z');

	cout << "Pressure X: " << pressureX << endl;
	cout << "Pressure Z: " << pressureZ << endl;

	//the spring value of this branch
	float spring = springConstant(b);

	//calculates the displacement value for each axis	
	float displacementX = pressureX / spring;
	float displacementZ = pressureZ / spring;

	//debug info
	cout << "length " << b->length << endl;
	cout << "Displacement - x: " << displacementX << "  z: " << displacementZ << endl;
	
	//mulitply a radians value by degrees variable to convert it from radians to degrees
	//not sure if this is needed...
	float degrees =  ( (float)(math::pi()) ) / 180.0f;

	//make sure no division of 0 is occuring 
	int len = b->length;
	if(len == 0){
		len = 0.00001f;
	}

	//currently returning NaN, asin() needs a value between [-1, 1] 
	//but currently getting a value too large
	float motionAngleX = asin(displacementX/float(len));
	float motionAngleZ = asin(displacementZ/float(len));

	cout << "Motion Angle - x: " << motionAngleX << "  z: " << motionAngleZ << endl;
	cout << endl;

	b->rotation.x = motionAngleX;
	b->rotation.z = motionAngleZ;

	//temporarily just rotating by displacement value because motionAngle is NaN
	b->rotation.x = ((displacementX*degrees) / 180 ) * 20;
	b->rotation.z = ((displacementZ*degrees) / 180 ) * 20;



	cout << "Restricted Angle - x: " << b->rotation.x << "  z: " << b->rotation.z << endl;
	cout << endl;

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

/*
	public method for toggling wind off and on
*/
void Tree::toggleWind(){
	windEnabled = ! windEnabled;
}

/* 
	sets the position we want to draw the tree at
*/
void Tree::setPosition(vec3 position) {
	m_position = position;
}

/* 
	sets the wind force
*/
void Tree::setWindForce(vec3 wind){
	windForce = wind;
}

void Tree::adjustWind(int axis, int dir){
	float increase = 1.0f;
	float aIncrease = 0.1f;

	if (axis == 'x'){
		if (dir == 1){
			windForce.x += increase;
		} else if (dir == -1){
			windForce.x -= increase;
		}
	} else if (axis == 'z'){
		if (dir == 1){
			windForce.z += increase;
		} else if (dir == -1){
			windForce.z -= increase;
		}
	}  else if (axis == 'a'){
		if (dir == 1){
			windCoefficent += aIncrease;
		} else if (dir == -1){
			windCoefficent -= aIncrease;
		}
	} 
}

/* Builds a test tree to work with for simulating wind animation.
	Trees is 'numBranches' segments tall, with 4 branches inbetween each segment.
*/
branch* Tree::makeDummyTree(int numBranches){
	//hardcoded values for this dummy tree
	float width = 0.3f;
	float length = 5.0f;

	branch* b = new branch();
	b->name = "trunk"+to_string(numBranches);
	b->direction = vec3(0,1,0);
	b->offset = math::random(0.0f,0.1f);
	b->length = length;
	b->baseWidth = width * numBranches;
	b->topWidth = (width * (numBranches - 1));
	if(numBranches == 1){
		b->topWidth = (width /2);
	}

	b->basisRot = vec3(0,0,0);
	if(numBranches > 1){

		for (int i = 0; i < 4; i++){
			branch* c = new branch();

			c->name = "branch"+to_string(i)+" trunk"+to_string(numBranches);

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