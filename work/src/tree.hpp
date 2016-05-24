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
	cgra::vec3 direction;
	cgra::vec3 basisRot;          // Euler angle rotations for the branch basis

	float length;

	float widthBase;
	float widthTop;

	cgra::vec3 rotation;          // Rotation of joint in the basis (degrees)

	std::vector<branch *> children;
};

class Tree{
	public:
		branch *root = nullptr;

		Tree();

		void renderTree();
		void setPosition(cgra::vec3);

	private:		
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);
		float rad = 0.015;

		branch* makeDummyTree(int);


		void renderBranch(branch *b);

		//drawing
		void drawBranch(branch*);
		void drawJoint();
		void drawAxis(branch*);
};