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

	std::string name;

	float length;
	float baseWidth;
	float topWidth;
	float offset;

	cgra::vec3 rotation = cgra::vec3(0,0,0);// Rotation of joint in the basis (degrees)

	std::vector<branch *> children = std::vector<branch*>();
	branch* parent;
};

class Tree{
	public:
		branch *root;

		Tree();

		void drawEnvelope();
		void renderTree();
		void renderStick();
		void renderAttractionPoints();
		
		void setPosition(cgra::vec3);
		void toggleWind();

	private:
		float param_branchLength = 0.4f;
		float param_radiusOfInfluence = 20 * param_branchLength;
		float param_killDistance = 2 * param_branchLength;

		float treeHeight;
		float trunkHeight;
		float maxX = 3.0f;
		float maxZ = 3.0f;
		float minX = -3.0f;
		float minZ = -3.0f;

		float yStep;
		float thetaStep = 15.0f;
		
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);
		cgra::vec3 windForce = cgra::vec3(0.0f, 0.0f, 0.0f);

		float width = 0.3f;
		float length = 5.0f;
		float elasticity = 20.0f;
		float time = 0.0f;
		bool windEnabled = false;

		std::vector<branch *> treeNodes;
		std::vector<std::vector<cgra::vec3>> envelope;
		std::vector<cgra::vec3> attractionPoints;

		branch* generateTree();
		float setWidth(branch*);
		std::vector<std::vector<int>> getAssociatedPoints();
		void cullAttractionPoints();
		void generateAttractionPoints(int num);
		void generateAttractionPointsVolumetric(int num);
		void generateEnvelope(int steps);
		float envelopeFunction(float u,float theta);

		bool inEnvelope(cgra::vec3);

		//drawing
		void renderBranch(branch *b);
		void drawBranch(branch*);
		void drawJoint(branch*);
		void drawAxis(branch*);

		void renderStick(branch *b);

		void setWindForce(cgra::vec3);
		float calculatePressure(branch*, float);
		float springConstant(branch*);
		void applyWind(branch*);
		float displacement(branch*, float);
};