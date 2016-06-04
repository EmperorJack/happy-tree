//-----------------------------
// 308 Final Project
// Class to represent the structure of a tree
//-----------------------------
#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

struct branch{
	//Position is at the end of the parent branch
	cgra::vec3 direction;			// initial direction of the branch
	cgra::vec3 basisRot;          	// Euler angle rotations for the branch basis

	std::string name;				// helpful with debug info

	float length;					// length of the branch
	float baseWidth;				// width of the base of the branch
	float topWidth;					// width of the top of the branch
	float offset;					// used to offset the branches sway in the wind

	cgra::vec3 rotation;          	// Rotation of joint in the basis (degrees)

	//branch* parent;
	std::vector<branch *> children;	// all child branches of this branch
};

class Tree{
	public:
		Tree();

		void renderTree();
		void setPosition(cgra::vec3);
		void toggleWind();

	private:		
		branch *root = nullptr; 	//the root section of the tree (first piece of trunk)

		//the position this tree will exist in world space
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);	

		branch* makeDummyTree(int);

		void renderBranch(branch *b);
		//drawing
		void drawBranch(branch*);
		void drawJoint(branch*);
		void drawAxis(branch*);

		// BELLOW HERE - wind variables and methods

		// the elasticity value of this tree, used for calculating spring value of branches
		float elasticity = 20.0f;	
		float time = 0.0f;			// time value used for moving along a sine curve 		
		bool windEnabled = false;
		//the wind acting upon this tree
		cgra::vec3 windForce = cgra::vec3(0.0f, 0.0f, 0.0f);	

		void setWindForce(cgra::vec3);
		float calculatePressure(branch*, float, int);
		float springConstant(branch*);
		void applyWind(branch*);
		float displacement(branch*, float);
};