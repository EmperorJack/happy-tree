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
	branch b;
	b.direction = vec3(0,1,0);
	b.length = 2.0f;
	b.widthBase = 0.5f;
	b.widthTop = 0.4f;
	root = &b;
}