#include "rendercontext3d.h"

#include <iostream>

extern short verbose;

// Render Context Creation Classes Implementation \\

// Render Context ABC Constructor
render_context::render_context(int sizex, int sizey)
	: wind_x(sizex), wind_y(sizey) {};

// Render Context OpenGL Constructor 
render_context_OGL::render_context_OGL(int xx, int yy, short v_maj, short v_min)
	: render_context(xx, yy), gl_ver_major(v_maj), gl_ver_minor(v_min)
{
	window_context();
	extensions_setup();
}

// Render Context OpenGL Destructor 
render_context_OGL::~render_context_OGL()
{
	glfwDestroyWindow(window);
	window = nullptr; 
	if (verbose) std::cout << "DBG::Main Render Window Destroyed Sucessfully \n \n";
}

// Render Context OpenGL - Window & Context Setup (GLFW). 
void render_context_OGL::window_context() 
{
	// GLFW Setup -
	glfwInit();
	if (!glfwInit())
	{
		std::cerr << "ERR::GLFW FAILED TO INITALIZE \n";
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, gl_ver_major);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, gl_ver_minor);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE); // Fixed Window Size. 
	glfwWindowHint(GLFW_SAMPLES, 4); // MSAA.

	window = glfwCreateWindow(wind_x, wind_y, "Stable Fluids 3D Solver | Render_API: OGL | Max_Threads: 8", NULL, NULL);

	if (window == NULL || window == nullptr)
	{
		std::cerr << "ERR::GLFW FAILED TO CREATE WINDOW \n";
		glfwTerminate();
	}

	glfwMakeContextCurrent(window); // In Current (main) thread. Where all Rendering Occurs. 
	glViewport(0, 0, wind_x, wind_y);

	std::cout << "DBG::GLFW Window and Initalzation Scuessful \n \n";
	if (verbose == 1)
	{
		std::cout << "DBG::GLFW Window and Initalzation Scuessful \n \n";
	}
}

// Render Context OpenGL - Extensions Setup (GLEW). 
void render_context_OGL::extensions_setup()
{
	// GLEW Setup
	glewExperimental = GL_TRUE;
	glewInit();
	if (glewInit() != GLEW_OK)
	{
		std::cerr << "ERR::GLEW FAILED TO INITALIZE \n";
	}

	// Query GL Device and Version Info - 
	render_device = glGetString(GL_RENDERER);
	version = glGetString(GL_VERSION);

	std::cout << "<OPENGL VERSION INFO BEGIN> \n";
	std::cout << "RENDER DEVICE = " << render_device << "\n";
	std::cout << "VERSION = " << version << "\n";
	std::cout << "<OPENGL VERSION INFO END> \n \n";
}

// Render Context OpenGL - Window Getter. 
GLFWwindow* render_context_OGL::get_window()
{
	return window; 
}



