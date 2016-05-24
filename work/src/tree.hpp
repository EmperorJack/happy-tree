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

		branch *root;
	private:
};