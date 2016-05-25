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
	float length;

	float widthBase;
	float widthTop;

	std::vector<branch *> children;
};

class Tree{
	public:
		Tree();
		void drawEnvelope();

		branch* root;
	private:
		std::vector<std::vector<cgra::vec3>> envelope;

		branch* makeDummyTree(int num);
		void generateEnvelope(int steps);
		float envelopeFunction(float u,float theta, float range);
};