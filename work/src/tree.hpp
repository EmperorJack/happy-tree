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
	cgra::vec3 position = cgra::vec3(0,0,0); //Only used while generating the tree
	cgra::vec3 direction = cgra::vec3(0,0,0);
	cgra::vec3 basisRot = cgra::vec3(0,0,0);// Euler angle rotations for the branch basis

	std::string name;				// helpful with debug info

	float maxX = 0.0f;
	float minX = 0.0f;
	float maxZ = 0.0f;
	float minZ = 0.0f;

	float length;					// length of the branch
	float baseWidth;				// width of the base of the branch
	float topWidth;					// width of the top of the branch
	float offset;					// used to offset the branches sway in the wind

	branch* parent;
	std::vector<branch *> children = std::vector<branch*>();	// all child branches of this branch

	cgra::vec3 worldDir = cgra::vec3(0,0,0);
	cgra::vec3 rotation = cgra::vec3(0,0,0);          // Rotation of joint in the basis (degrees)
};

class Tree{
	public:
		Tree();

		void drawEnvelope();
		void renderTree();
		void renderStick();
		void renderAttractionPoints();
		
		void setPosition(cgra::vec3);
		void toggleWind();
		void generateNewTree();
		void toggleTreeType();
		void adjustWind(int, int);	
		
	private:
		branch *root = nullptr; 	//the root section of the tree (first piece of trunk)
		branch *generatedTreeRoot = nullptr; 	//the root section of the tree (first piece of trunk)
		branch *dummyTreeRoot = nullptr; 	//the root section of the tree (first piece of trunk)

		//the position this tree will exist in world space
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);	

		float param_branchLength;
		float param_radiusOfInfluence;
		float param_killDistance;
		float param_branchTipWidth;
		float param_branchMinWidth;

		float treeHeight;
		float trunkHeight;
		float maxX = 3.0f;
		float maxZ = 3.0f;
		float minX = -3.0f;
		float minZ = -3.0f;

		float yStep;
		float thetaStep = 1.0f;

		std::vector<branch *> treeNodes;
		std::vector<std::vector<cgra::vec3>> envelope;
		std::vector<cgra::vec3> attractionPoints;

		//Tree generation Methods Start <<<<
		branch* generateTree();
		float setWidth(branch*);
		std::vector<std::vector<int>> getAssociatedPoints();
		void cullAttractionPoints();
		void generateAttractionPoints(int num);
		void generateAttractionPointsVolumetric(int num);
		void generateEnvelope(int steps);
		float envelopeFunction(float u,float theta);

		bool inEnvelope(cgra::vec3);
		branch* makeDummyTree(int);
		//Tree generation Methods STOP <<<<

		//drawing
		void renderBranch(branch *b, int depth=0);
		void drawBranch(branch*);
		void drawJoint(branch*);
		void drawAxis(branch*);
		void renderStick(branch *b, int depth=0);
		// BELLOW HERE - wind variables and methods

		// the elasticity value of this tree, used for calculating spring value of branches
		//hard coded right now, as number increases tree sway decreases
		//make skinnier trees have higher elasticity to prevent them from going crazy...
		float elasticity = 2000.0f;	
		float time = 0.0f;			// time value used for moving along a sine curve 		
		bool windEnabled = true;
		bool dummyTree = false;
		float windCoefficent = 1.2f;
		float timeIncrement = 0.004f;
		//the wind acting upon this tree
		cgra::vec3 desiredWindForce = cgra::vec3(0.0f, 0.0f, 0.0f);
		void setWindForce(cgra::vec3);
		
		float calculatePressure(branch*, float, int);
		float springConstant(branch*);
		void applyWind(branch*);

		void updateWorldWindDirection(branch*, cgra::vec3);
};