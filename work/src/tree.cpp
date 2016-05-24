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
	branch b;
	b.direction = vec3(0,1,0);
	b.length = 2.0f;
	b.widthBase = 0.5f * numBranches;
	b.widthTop = (0.5f * (numBranches - 1)) + 0.01f;
	if(numBranches > 1){
		b.children.push_back(makeDummyTree(numBranches - 1));
	}
	return &b;
}