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

struct treeNode{
	cgra::vec3 position;
	cgra::vec3 direction;
	float length;
	treeNode* parent;
};

class Tree{
	public:
		branch *root;

		Tree();

		void drawEnvelope();
		void renderTree();
		void renderAttractionPoints();
		
		void setPosition(cgra::vec3);

	private:
		float param_branchLength = 0.8f;
		float param_radiusOfInfluence = 17 * param_branchLength;
		float param_killDistance = 5 * param_branchLength;

		float treeHeight;
		float trunkHeight;
		float maxX = 3.0f;
		float maxZ = 3.0f;
		float minX = -3.0f;
		float minZ = -3.0f;

		float yStep;
		float thetaStep = 15.0f;
		
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);
		float width = 0.3f;
		float length = 5.0f;

		std::vector<treeNode *> treeNodes;
		std::vector<std::vector<cgra::vec3>> envelope;
		std::vector<cgra::vec3> attractionPoints;

		treeNode* generateTree();
		std::vector<std::vector<cgra::vec3>> getAssociatedPoints();
		void generateAttractionPoints(int num);
		void generateAttractionPointsVolumetric(int num);
		void generateEnvelope(int steps);
		float envelopeFunction(float u,float theta);

		bool inEnvelope(cgra::vec3);

		//drawing
		void renderBranch(branch *b);
		void drawBranch(branch*);
		void drawJoint();
		void drawAxis(branch*);

		branch* makeDummyTree(int);
};