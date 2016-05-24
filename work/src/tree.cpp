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
	root = makeDummyTree(4);
}

branch* Tree::makeDummyTree(int numBranches){
	branch* b = new branch();
	b->direction = vec3(0,1,0);
	b->length = 2.0f;
	b->widthBase = 0.5f * numBranches;
	b->widthTop = (0.5f * (numBranches - 1)) + 0.01f;
	b->basisRot = vec3(0,0,0);
	if(numBranches > 1){
		b->children.push_back(makeDummyTree(numBranches - 1));
	}
	return b;
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

			//draw the axes of this branch
			drawAxis(b);

			//draw the joint of this branch
			drawJoint();
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
		cgraSphere(1.2*rad);
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
		cgraCylinder(rad, rad/3, b->length);
	glPopMatrix();
}

void Tree::drawAxis(branch* b){
	//X-axes
	glPushMatrix();
		glColor3f(1,0,0);
		glRotatef(90,0, 1,0);
		cgraCylinder(0.3*rad, 0.3*rad, 4.0*rad);
		glTranslatef(0, 0, 4.0*rad);
		cgraCone(1.2*rad, 1.2*rad);
	glPopMatrix();

	//Y-axes
	glPushMatrix();
		glColor3f(0,1,0);
		glRotatef(-90,1,0,0);
		cgraCylinder(0.3*rad, 0.3*rad, 4.0*rad);
		glTranslatef(0,0, 4.0*rad);
		cgraCone(1.2*rad, 1.2*rad);
	glPopMatrix();

	//Z-axes
	glPushMatrix();
		glColor3f(0,0,1);
		cgraCylinder(0.3*rad, 0.3*rad, 4.0*rad);		
		glTranslatef(0,0,4.0*rad);
		cgraCone(1.2*rad, 1.2*rad);	
	glPopMatrix();
}

void Tree::setPosition(vec3 position) {
	m_position = position;
}