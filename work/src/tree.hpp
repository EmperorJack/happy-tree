//-----------------------------
// 308 Final Project
// Class to represent the structure of a tree
//-----------------------------
#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "fuzzy_object.hpp"

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

	Geometry* jointModel = nullptr;
	Geometry* branchModel = nullptr;
	FuzzyObject* branchFuzzySystem = nullptr;
};

class Tree{
	public:
		Tree(float height = 20.0f , float trunk = 0.0f, float branchLength = 2.0f ,float influenceRatio = 8.0f, float killRatio = 1.0f, float prm_branchTipWidth = 0.06f,float prm_branchMinWidth = 0.08f);
		~Tree();

		void drawEnvelope();
		void renderTree(bool);
		void renderStick();
		void renderAttractionPoints();

		void setPosition(cgra::vec3);
		void toggleWind();
		void toggleTreeType();
		void adjustWind(int, int);

		std::vector<Geometry*> getGeometries();

		// Fuzzy particle system methods
		void buildFuzzySystems(bool);
		bool finishedBuildingFuzzySystems();
		std::vector<cgra::vec3> getFuzzySystemPoints();

	private:
		branch* root = nullptr; 	//the root section of the tree (first piece of trunk)
		branch* generatedTreeRoot = nullptr; 	//the root section of the tree (first piece of trunk)
		branch* dummyTreeRoot = nullptr; 	//the root section of the tree (first piece of trunk)

		//the position this tree will exist in world space
		cgra::vec3 m_position = cgra::vec3(0.0f, 0.0f, 0.0f);

		float prm_branchLength;
		float prm_radiusOfInfluence;
		float prm_killDistance;
		float prm_branchTipWidth;
		float prm_branchMinWidth;

		float treeHeight;
		float trunkHeight;
		float maxX = 3.0f;
		float maxZ = 3.0f;
		float minX = -3.0f;
		float minZ = -3.0f;

		float yStep;
		float thetaStep = 10.0f;

		std::vector<branch *> treeNodes;
		std::vector<std::vector<cgra::vec3>> envelope;
		std::vector<cgra::vec3> attractionPoints;

		std::vector<FuzzyObject*> fuzzyBranchSystems;
		bool fuzzySystemFinishedBuilding = false;

		//Tree Generation Methods
		branch* generateTree();
		float setWidth(branch*);
		void generateGeometry(branch*);
		std::vector<std::vector<int>> getAssociatedPoints();
		void cullAttractionPoints();
		void generateAttractionPoints(int num);
		void generateAttractionPointsVolumetric(int num);
		void generateEnvelope(int steps);
		float envelopeFunction(float u,float theta);

		void getBranchGeometry(branch*, std::vector<Geometry*>*);
		void getBranchFuzzySystemPoints(branch*, std::vector<cgra::vec3>*);

		bool inEnvelope(cgra::vec3);
		branch* makeDummyTree(int);

		//Drawing Methods
		void renderBranch(branch *b, bool, int depth=0);
		void drawBranch(branch*, bool);
		void drawJoint(branch*, bool);
		void drawAxis(branch*);
		void renderStick(branch *b, int depth=0);

		//Wind Simulation

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

		cgra::mat3 angleAxisRotation(float, cgra::vec3);
};
