//---------------------------------------------------------------------------
// COMP308 Final Project
//---------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>

#include "cgra_geometry.hpp"
#include "cgra_math.hpp"
#include "simple_image.hpp"
#include "simple_shader.hpp"
#include "simple_gui.hpp"
#include "opengl.hpp"
#include "geometry.hpp"
#include "tree.hpp"
#include "fuzzy_object.hpp"
#include "particle_system.hpp"

using namespace std;
using namespace cgra;

// Window
GLFWwindow* g_window;

// Frame related values
int frameCount = 0;
double frameRate = 0.0;

// Projection values
float g_fovy = 20.0;
float g_znear = 0.1;
float g_zfar = 1000.0;

// Mouse controlled values
bool g_leftMouseDown = false;
bool g_rightMouseDown = false;
vec2 g_mousePosition;
float g_pitch = 0;
float g_yaw = 0;
float g_zoom = 1.0;

// Shader fields
GLuint g_shader = 0;

// Geometry draw lists
Geometry* g_model = nullptr;
Geometry* g_terrain = nullptr;

// Tree to animate
Tree* g_tree = nullptr;
float tree_h = 20.0f;
float tree_t = 1.0f;
float tree_bL = 2.0f;
float tree_inf = 8.0f;
float tree_kill = 1.0f;
float tree_tW = 0.04;
float tree_mW = 0.08;

// Particle system fields
FuzzyObject* g_fuzzy_system = nullptr;
ParticleSystem* g_treeParticleSystem = nullptr;
bool treeFuzzySystemFinishedBuilding = false;

float spawnPointShiftAmount = 0.1f;
bool explodingSystem = false;

// Toggle fields
bool drawAxes = true;
bool treeMode = false;
bool wireframeMode = false;
bool realtimeBuild = false;
bool partyMode = false;

// Mouse Position callback
void cursorPosCallback(GLFWwindow* win, double xpos, double ypos) {
	// cout << "Mouse Movement Callback :: xpos=" << xpos << "ypos=" << ypos << endl;
	if (g_leftMouseDown) {
		g_yaw -= g_mousePosition.x - xpos;
		g_pitch -= g_mousePosition.y - ypos;
	}
	g_mousePosition = vec2(xpos, ypos);
}

// Mouse Button callback
void mouseButtonCallback(GLFWwindow *win, int button, int action, int mods) {
	// cout << "Mouse Button Callback :: button=" << button << "action=" << action << "mods=" << mods << endl;
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		g_leftMouseDown = (action == GLFW_PRESS);
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		g_rightMouseDown = (action == GLFW_PRESS);
		if (g_rightMouseDown) {
			g_tree->buildFuzzySystems(true);
		}
	}

	if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
		if (action == GLFW_PRESS) {
			for (int i = 0; i < 100; i++) {
				g_tree->buildFuzzySystems(true);
			}
		}
	}
}

// Scroll callback
void scrollCallback(GLFWwindow *win, double xoffset, double yoffset) {
	// cout << "Scroll Callback :: xoffset=" << xoffset << "yoffset=" << yoffset << endl;
	g_zoom -= yoffset * g_zoom * 0.2;
}

void initLight();

// Keyboard callback
void keyCallback(GLFWwindow *win, int key, int scancode, int action, int mods) {
	//cout << "Key Callback :: key=" << key << "scancode=" << scancode	<< "action=" << action << "mods=" << mods << endl;
	//Tree Gen Stuff
	if (mods == 2) {
		if (key == 'R' && action == 1) {
			delete(g_tree);
			g_tree = new Tree(tree_h, tree_t, tree_bL, tree_inf, tree_kill, tree_tW, tree_mW);

			delete(g_treeParticleSystem);
			treeFuzzySystemFinishedBuilding = false;
			realtimeBuild = false;
		}
		if (key == 'T' && action == 1) {
			treeMode = !treeMode;
		}

		if (key == 'W' && action == 1) {
			tree_h += 1.0f;
		}
		if (key == 'S' && action == 1) {
			tree_h -= 1.0f;
			tree_h = tree_h < 1.0 ? 1.0 : tree_h;
		}

		if (key == 'D' && action == 1) {
			tree_bL += 0.1f;
		}
		if (key == 'A' && action == 1) {
			tree_bL -= 0.1f;
			tree_bL = tree_bL < 1.0 ? 1.0 : tree_bL;
		}

		if (key == 'E' && action == 1) {
			tree_mW += 0.01f;
		}
		if (key == 'Q' && action == 1) {
			tree_mW -= 0.01f;
			tree_mW = tree_mW <= 0 ? 0.01 : tree_mW;
		}
	}else{
		if (key == 'A' && action == 1) {
			drawAxes = !drawAxes;
		}

		if (key == 'F' && action == 1) {
			g_tree->toggleWind();
		}

		if (key == 'C' && action == 1) {
			g_tree->toggleTreeType();
		}

		// increase wind on X axis
		if (key == 'J' && (action == 1 || action == 2)) {
			g_tree->adjustWind('x', 1);
		}

		// decrease wind on X axis
		if (key == 'N' && (action == 1 || action == 2)) {
			g_tree->adjustWind('x', -1);
		}

		// increase wind on Z axis
		if (key == 'K' && (action == 1 || action == 2)) {
			g_tree->adjustWind('z', 1);
		}

		// decrease wind on Z axis
		if (key == 'M' && (action == 1 || action == 2)) {
			g_tree->adjustWind('z', -1);
		}

		// increase a coefficient in wind calculation
		if (key == 'H' && (action == 1 || action == 2)) {
			g_tree->adjustWind('a', 1);
		}

		// decrease a coefficient in wind calculation
		if (key == 'B' && (action == 1 || action == 2)) {
			g_tree->adjustWind('a', -1);
		}

		// increase a coefficient in wind calculation
		if (key == 'G' && (action == 1 || action == 2)) {
			g_tree->adjustWind('t', 1);
		}

		// decrease a coefficient in wind calculation
		if (key == 'V' && (action == 1 || action == 2)) {
			g_tree->adjustWind('t', -1);
		}

		// 'p' key pressed
		if (key == 'P' && action == 1) {
			partyMode = !partyMode;

			// Reset the light properties
			initLight();
		}

		// 'w' key pressed
		if (key == 'W' && action == 1) {
			wireframeMode = !wireframeMode;
		}

		// 'space' key pressed
		if (key == 32 && action == 1) {
			explodingSystem = !explodingSystem;

			// Check if the tree finished building it's particle systems
			if (g_tree->finishedBuildingFuzzySystems()) {
				g_treeParticleSystem = new ParticleSystem(g_tree->getFuzzySystemPoints());
				treeFuzzySystemFinishedBuilding = true;
			 	g_treeParticleSystem->explode();
			}
		}

		// 'r' key pressed
		if (key == 'R' && (action == 1)) {
			explodingSystem = false;
			if (treeFuzzySystemFinishedBuilding) g_treeParticleSystem->resetParticles();
		}

		// 'q' key pressed
		if (key == 'Q' && action == 1) {
			realtimeBuild = !realtimeBuild;
		}

		// 'v' key pressed
		if (key == 'L' && action == 1) {
			g_fuzzy_system->toggleParticleViewMode();
		}
	}
}

// Character callback
void charCallback(GLFWwindow *win, unsigned int c) {
	// cout << "Char Callback :: c=" << char(c) << endl;
}

// Load and setup the 3D geometry models
void initGeometry() {
	//g_model = new Geometry("./work/res/assets/sphere.obj");
	//g_model = generateSphereGeometry(3.0f, 6, 6);
	g_model = generateCylinderGeometry(1.0f, 1.0f, 5.0f, 4, 4);
	g_model->setPosition(vec3(5, 1, 5));

	g_terrain = new Geometry("./work/res/assets/plane.obj");

	g_tree = new Tree();
	g_tree->setPosition(vec3(0, 0, 0));
}

// Setup the materials per geometric object
void initMaterials() {
	vec4 black = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 grey = vec4(0.2, 0.2, 0.2, 1.0);
	vec4 white = vec4(1.0, 1.0, 1.0, 1.0);
	
	g_terrain->setMaterial(vec4(0.5,0.5,0.5,1.0), vec4(0.5,0.5,0.5,1.0), vec4(0.1,0.1,0.1,1.0), 20.0f,black);

	g_model->setMaterial(grey, vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, black);
}

// Loads in a texture from the given location
GLuint initTexture(string path) {
	GLuint g_texture;
	Image tex(path);

	glActiveTexture(GL_TEXTURE0); // Use slot 0, need to use GL_TEXTURE1 ... etc if using more than one texture PER OBJECT
	glGenTextures(1, &g_texture); // Generate texture ID
	glBindTexture(GL_TEXTURE_2D, g_texture); // Bind it as a 2D texture

	// Setup sampling strategies
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// Enable repeating of the texture
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// Finnaly, actually fill the data into our texture
	gluBuild2DMipmaps(GL_TEXTURE_2D, 3, tex.w, tex.h, tex.glFormat(), GL_UNSIGNED_BYTE, tex.dataPointer());

	return g_texture;
}

// Sets up the lights and their individual properties
void initLight() {
	vec4 black = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 white = vec4(1.0, 1.0, 1.0, 1.0);

	// Point light
	glLightfv(GL_LIGHT0, GL_AMBIENT, black.dataPointer());
	glLightfv(GL_LIGHT0, GL_DIFFUSE, vec4(0.75, 0.75, 1.0, 1.0).dataPointer());
	glLightfv(GL_LIGHT0, GL_SPECULAR, white.dataPointer());

	// Point light
	glLightfv(GL_LIGHT1, GL_AMBIENT, black.dataPointer());
	glLightfv(GL_LIGHT1, GL_DIFFUSE, vec4(1.0, 0.75, 0.75, 1.0).dataPointer());
	glLightfv(GL_LIGHT1, GL_SPECULAR, white.dataPointer());
}

// Load a shader from a given location
void initShader(string vertPath, string fragPath) {
	g_shader = makeShaderProgramFromFile({GL_VERTEX_SHADER, GL_FRAGMENT_SHADER }, { vertPath, fragPath });
}

// Sets up where the camera is in the scene
void setupCamera(int width, int height) {
	// Set up the projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(g_fovy, width / float(height), g_znear, g_zfar);

	// Set up the view part of the model view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0, 0, -50 * g_zoom);
	glRotatef(g_pitch, 1, 0, 0);
	glRotatef(g_yaw, 0, 1, 0);
}

// Sets up the lighting of the scene
void setupLight() {
	// Set the light positions and directions
	glLightfv(GL_LIGHT0, GL_POSITION, vec4(5.0, 5.0, 5.0, 1.0).dataPointer());
	glLightfv(GL_LIGHT1, GL_POSITION, vec4(-5.0, 5.0, -5.0, 1.0).dataPointer());

	if (frameCount % 20 == 0 && partyMode) {
		glLightfv(GL_LIGHT0, GL_DIFFUSE, vec4(math::random(0.0f, 1.0f), math::random(0.0f, 1.0f), math::random(0.0f, 1.0f), 1.0).dataPointer());
		glLightfv(GL_LIGHT1, GL_DIFFUSE, vec4(math::random(0.0f, 1.0f), math::random(0.0f, 1.0f), math::random(0.0f, 1.0f), 1.0).dataPointer());
	}
}

// Render the global axes for reference
void renderGlobalAxes() {
	glDisable(GL_LIGHTING);
	glLineWidth(2);

	// X axis
	glColor3f(1.0f, 0.5f, 0.5f);
	glBegin(GL_LINES);
	glVertex3f(-100, 0, 0);
	glVertex3f(100, 0, 0);
	glEnd();

	// Y axis
	glColor3f(0.5f, 1.0f, 0.5f);
	glBegin(GL_LINES);
	glVertex3f(0, -100, 0);
	glVertex3f(0, 100, 0);
	glEnd();

	// Z axis
	glColor3f(0.5f, 0.5f, 1.0f);
	glBegin(GL_LINES);
	glVertex3f(0, 0, -100);
	glVertex3f(0, 0, 100);
	glEnd();

	glEnable(GL_LIGHTING);
}

// Render a square plane with the given length
void renderPlane(float length) {
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vec4(0.75f, 0.75f, 0.75f, 1.0f).dataPointer());
	glMaterialfv(GL_FRONT, GL_SPECULAR, vec4(0.5f, 0.5f, 0.5f, 1.0f).dataPointer());
	float shininess = 128.0f;
	glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	float l = length / 2.0f;

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glBegin(GL_QUADS);
		glNormal3f(0, 1, 0);
		glTexCoord2f(-1, -1);
		glVertex3f(-l, 0, -l);
		glTexCoord2f(-1, 1);
		glVertex3f(-l, 0, l);
		glTexCoord2f(1, 1);
		glVertex3f(l, 0, l);
		glTexCoord2f(1, -1);
		glVertex3f(l, 0, -l);
	glEnd();
}

void update() {
	if (!treeFuzzySystemFinishedBuilding) {

		// Update tree particle system building
		if (!g_tree->finishedBuildingFuzzySystems() && realtimeBuild) {
			g_tree->buildFuzzySystems(true);
		}
	} else {
		g_treeParticleSystem->update();
	}
}

// Render the geometry of the scene
void renderScene() {
	if (partyMode) glRotatef(frameCount * -1.5f, 0, 1, 0);

	// Render terrain
	g_terrain->renderGeometry(false);

	if (treeMode){
		glDisable(GL_LIGHTING);
		//Render Tree
		g_tree->renderStick();
		g_tree->drawEnvelope();
		// g_tree->renderAttractionPoints();
		glEnable(GL_LIGHTING);
	} else {

		// Render geometry
		glUniform1i(glGetUniformLocation(g_shader, "useTexture"), true);
		g_tree->renderTree(wireframeMode);
		glUniform1i(glGetUniformLocation(g_shader, "useTexture"), false);
	}

	// Render complete fuzzy particle system
	if (treeFuzzySystemFinishedBuilding) {
		g_treeParticleSystem->render();
	}
}

// Draw the scene
void render(int width, int height) {

	// Update mechanics in the scene
	update();

	glViewport(0, 0, width, height);

	// Clear the background
	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	if (partyMode) glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	// Enable flags for normal rendering
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	// Render global axes
	if (drawAxes) renderGlobalAxes();

	// Setup the lighting, camera and shader
	setupLight();
	setupCamera(width, height);
	glUseProgram(g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "texture0"), 0);

	renderScene();

	// Disable flags for cleanup
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glUseProgram(0);
}

// Render the IMGUI overlay
void renderGUI() {
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// Start registering GUI components
	SimpleGUI::newFrame();

	// FPS counter
	char fpsString[128];
	sprintf(fpsString, "FPS: %.2f", 1 / frameRate);

	ImGui::SetNextWindowPos(ImVec2(10, 10));
	ImGui::Begin("", nullptr, ImVec2(0, 0), 0.3f,
				 ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoSavedSettings);

	ImGui::Text(string(fpsString).c_str());
	ImGui::Text(("Particle Count: " + to_string(g_tree->getFuzzySystemParticleCount())).c_str());

	ImGui::End();

	// Flush components and render
	SimpleGUI::render();
}

void APIENTRY debugCallbackARB(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar*, GLvoid*);

// Main program
int main(int argc, char **argv) {

	// Initialize the GLFW library
	if (!glfwInit()) {
		cerr << "Error: Could not initialize GLFW" << endl;
		abort();
	}
	int glfwMajor, glfwMinor, glfwRevision;
	glfwGetVersion(&glfwMajor, &glfwMinor, &glfwRevision);

	// Create a windowed mode window and its OpenGL context
	g_window = glfwCreateWindow(1024, 768, "Project", nullptr, nullptr);
	if (!g_window) {
		cerr << "Error: Could not create GLFW window" << endl;
		abort();
	}
	glfwMakeContextCurrent(g_window);

	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);

	// Initialize GLEW
	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		cerr << "Error: " << glewGetErrorString(err) << endl;
		abort();
	}

	// Print out our OpenGL verisions
	cout << "Using OpenGL " << glGetString(GL_VERSION) << endl;
	cout << "Using GLEW " << glewGetString(GLEW_VERSION) << endl;
	cout << "Using GLFW " << glfwMajor << "." << glfwMinor << "." << glfwRevision << endl;

	// Attach input callbacks to g_window
	glfwSetCursorPosCallback(g_window, cursorPosCallback);
	glfwSetMouseButtonCallback(g_window, mouseButtonCallback);
	glfwSetScrollCallback(g_window, scrollCallback);
	glfwSetKeyCallback(g_window, keyCallback);
	glfwSetCharCallback(g_window, charCallback);

	// Enable GL_ARB_debug_output if available. Not nessesary, just helpful
	if (glfwExtensionSupported("GL_ARB_debug_output")) {
		// This allows the error location to be determined from a stacktrace
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		// Set the up callback
		glDebugMessageCallbackARB(debugCallbackARB, nullptr);
		glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, true);
		cout << "GL_ARB_debug_output callback installed" << endl;
	}
	else {
		cout << "GL_ARB_debug_output not available. No worries." << endl;
	}

	// Initialize IMGUI
	if (!SimpleGUI::init(g_window, false)) {
		cerr << "Error: Could not initialize IMGUI" << endl;
		abort();
	}

	// Initialize geometry, materials, lighting and shaders
	initGeometry();
	initMaterials();
	initLight();
	initShader("./work/res/shaders/phongShader.vert", "./work/res/shaders/phongShader.frag");
	initTexture("./work/res/textures/bark.png");

	g_fuzzy_system = new FuzzyObject(g_model);

	double lastTime = glfwGetTime();
	int framesThisSecond = 0;

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(g_window)) {

		// FPS update
		double currentTime = glfwGetTime();
		framesThisSecond++;
		if (currentTime - lastTime >= 1.0) {
	    	frameRate = 1.0 / double(framesThisSecond);
	    	framesThisSecond = 0;
	    	lastTime += 1.0;
	    }

		int width, height;
		glfwGetFramebufferSize(g_window, &width, &height);

		// Main render
		render(width, height);

		// Render GUI on top
		renderGUI();

		glfwSwapBuffers(g_window);
		glfwPollEvents();

		// Update total frame count
		frameCount++;
	}

	glfwTerminate();
}

//-------------------------------------------------------------
// Fancy debug stuff
//-------------------------------------------------------------

// Function to translate source to string
string getStringForSource(GLenum source) {
	switch (source) {
	case GL_DEBUG_SOURCE_API:
		return("API");
	case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
		return("Window System");
	case GL_DEBUG_SOURCE_SHADER_COMPILER:
		return("Shader Compiler");
	case GL_DEBUG_SOURCE_THIRD_PARTY:
		return("Third Party");
	case GL_DEBUG_SOURCE_APPLICATION:
		return("Application");
	case GL_DEBUG_SOURCE_OTHER:
		return("Other");
	default:
		return("n/a");
	}
}

// Function to translate severity to string
string getStringForSeverity(GLenum severity) {
	switch (severity) {
	case GL_DEBUG_SEVERITY_HIGH:
		return("HIGH!");
	case GL_DEBUG_SEVERITY_MEDIUM:
		return("Medium");
	case GL_DEBUG_SEVERITY_LOW:
		return("Low");
	default:
		return("n/a");
	}
}

// Function to translate type to string
string getStringForType(GLenum type) {
	switch (type) {
	case GL_DEBUG_TYPE_ERROR:
		return("Error");
	case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
		return("Deprecated Behaviour");
	case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
		return("Undefined Behaviour");
	case GL_DEBUG_TYPE_PORTABILITY:
		return("Portability Issue");
	case GL_DEBUG_TYPE_PERFORMANCE:
		return("Performance Issue");
	case GL_DEBUG_TYPE_OTHER:
		return("Other");
	default:
		return("n/a");
	}
}

void APIENTRY debugCallbackARB(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei, const GLchar* message, GLvoid*) {
	if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;

	cerr << endl;
	cerr << "Type: " <<
		getStringForType(type) << "; Source: " <<
		getStringForSource(source) << "; ID: " << id << "; Severity: " <<
		getStringForSeverity(severity) << endl;
	cerr << message << endl;

	if (type == GL_DEBUG_TYPE_ERROR_ARB) throw runtime_error("");
}
