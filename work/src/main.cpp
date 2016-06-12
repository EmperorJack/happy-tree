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
Geometry* sphere = nullptr;
Geometry* cylinder = nullptr;

// Tree to animate
Tree* g_tree = nullptr;

// Particle system fields
FuzzyObject* g_fuzzy_system = nullptr;
float spawnPointShiftAmount = 0.1f;
bool explodingSystem = false;

// Toggle fields
bool drawAxes = false;
bool treeMode = false;
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
			g_fuzzy_system->buildSystemIncrement();
		}
	}

	if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
		if (action == GLFW_PRESS) {
			for (int i = 0; i < 100; i++) g_fuzzy_system->buildSystemIncrement();
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
	// cout << "Key Callback :: key=" << key << "scancode=" << scancode
	// 	<< "action=" << action << "mods=" << mods << endl;

	// 'a' key pressed
	if (key == 'A' && action == 1) {
		drawAxes = !drawAxes;
	}

	if (key == 'T' && action == 1) {
		treeMode = !treeMode;
	}

	// 'p' key pressed
	if (key == 'P' && action == 1) {
		partyMode = !partyMode;

		// Reset the light properties
		initLight();
	}

	// 'w' key pressed
	if (key == 'W' && action == 1) {
		g_model->toggleWireframe();
		cylinder->toggleWireframe();
		sphere->toggleWireframe();
	}

	// 'space' key pressed
	if (key == 32 && action == 1) {
		if (g_fuzzy_system->finishedBuilding()) {
			g_fuzzy_system->explode();
			explodingSystem = true;
		}
	}

	// 'up' key pressed
	if (key == 265 && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.z += spawnPointShiftAmount;
	}

	// 'left' key pressed
	if (key == 263 && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.x += spawnPointShiftAmount;
	}

	// 'right' key pressed
	if (key == 262 && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.x -= spawnPointShiftAmount;
	}

	// 'down' key pressed
	if (key == 264 && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.z -= spawnPointShiftAmount;
	}

	// 'e' key pressed
	if (key == 'E' && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.y -= spawnPointShiftAmount;
	}

	// 'r' key pressed
	if (key == 'R' && (action == 1 || action == 2)) {
		g_fuzzy_system->spawnPoint.y += spawnPointShiftAmount;
	}

	// 'q' key pressed
	if (key == 'Q' && action == 1) {
		realtimeBuild = !realtimeBuild;
	}

	// 'v' key pressed
	if (key == 'V' && action == 1) {
		g_fuzzy_system->toggleParticleViewMode();
	}
}

// Character callback
void charCallback(GLFWwindow *win, unsigned int c) {
	// cout << "Char Callback :: c=" << char(c) << endl;
}

// Load and setup the 3D geometry models
void initGeometry() {
	g_model = new Geometry("./work/res/assets/ignored/bunny-reduced.obj");
	g_model->setPosition(vec3(0, 0, 0));

	g_tree = new Tree();
	g_tree->setPosition(vec3(0, 0, 0));

	cylinder = generateCylinderGeometry(1.0f, 1.0f, 5.0f, 4, 4);
	cylinder->setPosition(vec3(0, 3, 0));

	sphere = generateSphereGeometry(2.0f, 4, 4);
	sphere->setPosition(vec3(0, 3, 0));
}

// Setup the materials per geometric object
void initMaterials() {
	vec4 black = vec4(0.0, 0.0, 0.0, 1.0);
	vec4 grey = vec4(0.2, 0.2, 0.2, 1.0);
	vec4 white = vec4(1.0, 1.0, 1.0, 1.0);

	g_model->setMaterial(grey, vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, black);

	cylinder->setMaterial(grey, vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, black);
	sphere->setMaterial(grey, vec4(0.8, 0.8, 0.8, 1.0), vec4(0.8, 0.8, 0.8, 1.0), 128.0f, black);
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

// Render the geometry of the scene
void renderScene() {
	if (partyMode) glRotatef(frameCount * -1.5f, 0, 1, 0);

	// Render plane
	//renderPlane(20);

	if (treeMode){
		//Render Tree
		glDisable(GL_LIGHTING);
		g_tree->renderTree();
		glEnable(GL_LIGHTING);
	} else {
		// Render geometry
		if (!g_fuzzy_system->finishedBuilding()) g_model->renderGeometry();
	}

	// Update building particle system
	if (realtimeBuild && !g_fuzzy_system->finishedBuilding()) g_fuzzy_system->buildSystemIncrement();

	if (explodingSystem && g_fuzzy_system->finishedBuilding()) g_fuzzy_system->updateSystem();

	// Render particle system
	g_fuzzy_system->renderSystem();

	//cylinder->renderGeometry();
	//sphere->renderGeometry();
}

// Draw the scene
void render(int width, int height) {
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
	ImGui::Text(("Particle Count: " + to_string(g_fuzzy_system->getParticleCount())).c_str());

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
