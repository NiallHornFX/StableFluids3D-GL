#include "renderobject3d.h"
#include "fluidobj3d.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

extern short verbose;

// Render Object Creation Classes Implementation - 

// Oppose to 2D Window does NOT Correspond to Grid in Size or Cell-Pixel Mapping. We are defining a window of aribitrary size, whose fragments will raymarch
// the FluidObject Grids passed to GPU as 3D Textures sampled in WS. 

// RenderObject_3D - ABC Implementations - 

// Render Object 3D ABC Constructor
renderobject_3D::renderobject_3D(const char *api_name, int v_maj, int v_min, const vec2<int> &win_size, const vec3<int> &grd_size, short rmode)
	: API_Name(api_name), ver_major(v_maj), ver_minor(v_min), window_size(win_size), grid_size(grd_size), rendermode(rmode) {}

// RenderObject_2D_OGL (OpenGL API) Implementations - 

// RenderObject_2D_OGL Constructor - 
renderobject_3D_OGL::renderobject_3D_OGL(const char *api_name, int v_maj, int v_min, const vec2<int> &ws, const vec3<int> &gs, GLFWwindow *winptr, short rmode)
	: renderobject_3D(api_name, v_maj, v_min, ws, gs, rmode), window_ptr(winptr)  // Initalize ABC Members Via Its Own Constructor. 
{
	std::cout << "DBG::RenderObject_3D Created For Render API: " << api_name << " " << v_maj << "." << v_min << "\n \n";
	// Call Setup MFuncs

	int v_ret = vertex_setup(); 
	int s_ret = shader_loader("shaders/render3d_ogl_vertShader.vert", "shaders/render3d_ogl_fragShader.frag" ); // HC ShaderPaths. 

	// Gen 3D Textures -
	glGenTextures(1, &tex_dens);
	glGenTextures(1, &tex_vel);
	// Do Not Bind On Initalization.
}

// RenderObject_2D_OGL Destructor -
renderobject_3D_OGL::~renderobject_3D_OGL()
{
	if (vert_shader || vert_shader != NULL) {
		glDeleteShader(vert_shader);
		vert_shader = NULL;
	}

	if (frag_shader || frag_shader != NULL) {
		glDeleteShader(frag_shader);
		frag_shader = NULL;
	}

	if (shader_prog || shader_prog != NULL) {
		glDeleteProgram(shader_prog);
		shader_prog = NULL; 
	}

	delete vertices; vertices = nullptr; 
	delete indices; indices = nullptr; 

	delete vert_shader_code; vert_shader_code = nullptr; 
	delete frag_shader_code; frag_shader_code = nullptr; 
}

// RenderObject_2D_OGL Vertex Setup Implementation, For Quad and Indexed Vertex Drawing - 
int renderobject_3D_OGL::vertex_setup()
{
	// Vertex Screen Quad Setup - 

	vertices = new GLfloat[12]
	{
		1.0f, 1.0f, 0.0f,    // 0 TR Vert
		1.0f, -1.0f, 0.0f,   // 1  BR Vert
	    -1.0f, -1.0f, 0.0f,  // 2 BL Vert
		-1.0f, 1.0f, 0.0f    // 3 TL Vert
	};

	indices = new GLuint[6]
	{
		0,1,3, // Tri 1
		1,2,3  // Tri 2
	};

	// Vertex Array/Buffer Setup Gen/Bind/AttribData - 
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 12, vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0); // VertPos
	glEnableVertexAttribArray(0);

	// Element Index Buffer Setup - 
	glGenBuffers(1, &EBO);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 6, indices, GL_STATIC_DRAW);

	return 0;
}

// RenderObject_3D_OGL Shader Compilation Checker MFunc Implementation - 
void renderobject_3D_OGL::shader_checkCompile(const char *type)
{
	assert(type == (const char*) "vertex" || type == (const char*) "fragment");

	const int len = 512; // length of error log char array. 
	int sucess_v, sucess_f;
	char err_log_v[len], err_log_f[len];

	// Vertex Shader Check - 
	if (strcmp(type, "vertex") == 0)
	{
		glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &sucess_v);
		if (!sucess_v)
		{
			glGetShaderInfoLog(vert_shader, len, NULL, err_log_v);
			std::cerr << "ERR::VERTEX SHADER COMPILE FAILED \n " << err_log_v << std::endl; 
		}
		else if (verbose == 1)
		{
			std::cout << "DBG::VERTEX SHADER COMPILE SUCEEDED \n";
		}
	}
	
	// Fragment Shader Check - 
	if (strcmp(type, "fragment") == 0)
	{
		glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &sucess_f);
		if (!sucess_f)
		{
			glGetShaderInfoLog(frag_shader, len, NULL, err_log_f);
			std::cerr << "ERR::FRAGMENT SHADER COMPILE FAILED \n " << err_log_f << std::endl;
		}
		else if (verbose == 1)
		{
			std::cout << "DBG::FRAGMENT SHADER COMPILE SUCEEDED \n";
		}
	}
}

// RenderObject_3D_OGL ShaderProgram Linker Checker MFunc Implementation - 
void renderobject_3D_OGL::shader_checkLink()
{
	int sucess;
	const int len = 512;
	char err_log[len];

	glGetProgramiv(shader_prog, GL_LINK_STATUS, &sucess);

	// If Not Sucess (thus sucess undefined/non-intialzed) print error log out (via cout)
	if (!sucess)
	{
		glGetProgramInfoLog(shader_prog, len, NULL, err_log);
		std::cout << "ERR:SHADER-PROGRAM: " << "LINKAGE_FAILED" << std::endl;
		std::cout << err_log << std::endl;
	}
	else if (verbose == 1)
	{
		std::cout << "DBG::SHADER:PROGRAM " << shader_prog << " LINKAGE SUCEEDED" << std::endl;
	}
}

// RenderObject_3D_OGL Shader Loader Implementation - 
int renderobject_3D_OGL::shader_loader(const char *vert_path, const char *frag_path)
{
	std::ifstream vert_shader_load, frag_shader_load; 
	std::stringstream v_shad_buf, f_shad_buf; 
	std::string temp_v, temp_f; 

	// Set ifstream exceptions
	vert_shader_load.exceptions(std::ios::badbit | std::ios::failbit);
	frag_shader_load.exceptions(std::ios::badbit | std::ios::failbit);

	// Shader Load And Write to Code Buffer. (Probs could of done this with strcpy directly much easier !) 
	// No Parsing, (For now) just check compile and linkage on shader creation below. 
	try
	{
		vert_shader_load.open(vert_path);
		frag_shader_load.open(frag_path);

		v_shad_buf << vert_shader_load.rdbuf(); temp_v = v_shad_buf.str();
		f_shad_buf << frag_shader_load.rdbuf(); temp_f = f_shad_buf.str();

		vert_shader_load.close(); 
		frag_shader_load.close(); 

		if (vert_shader_load.is_open() || frag_shader_load.is_open()) std::cerr << "ERR::Shader Closed Incorrectly \n";
	}
	catch (std::ifstream::failure err)
	{
		std::cerr << "ERR::Shader Load Err: " << err.what() << "\n";
		std::abort();
		return 1; 
	}

	// Fill Shader Code Member CStrs. 
	vert_shader_code = temp_v.c_str(); 
	frag_shader_code = temp_f.c_str(); 

	// Shader Compilation & Link - 
	vert_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader, 1, &vert_shader_code, NULL);
	glCompileShader(vert_shader);

	frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader, 1, &frag_shader_code, NULL);
	glCompileShader(frag_shader);

	shader_checkCompile("vertex"); shader_checkCompile("fragment");

	// Shader Program - 
	shader_prog = glCreateProgram();
	glAttachShader(shader_prog, vert_shader); glAttachShader(shader_prog, frag_shader);
	glLinkProgram(shader_prog);

	shader_checkLink();

	return 0; 
}

void renderobject_3D_OGL::shader_pipe(fluidobj_3d *f3obj)
{
	// UNIFORM CONSTANTS (Per Step) \\

	// Pass Window Size to GL Uniform. // Could be done after Shader Compile (Once.?)
	glUniform1i(glGetUniformLocation(shader_prog, "W_Size"), window_size.x); // Window Size (Assume Square Dim, Used for Frag-UV Space)

	// Interactive Render Mode Switching - 
	int ent_state = glfwGetKey(window_ptr, GLFW_KEY_ENTER);
	// Switch Between 0 and 1 (Density or Vel Render).
	if (ent_state == GLFW_PRESS)
	{
		// Switch Render Mode On Enter Press (Assuming 0 or 1 Modes Only).
		rendermode = !rendermode;

		// Update RenderMode Uniform - 
		glUniform1i(glGetUniformLocation(shader_prog, "Mode"), rendermode);
		if (verbose) std::cout << "DBG::RENDER MODE SWITCHED = " << rendermode << "\n";
	}

	// Current Step Uniform (RenderObj BaseMember et passed in solvestep). Fix naughty Double-Int cast. 
	glUniform1i(glGetUniformLocation(shader_prog, "Step"), (GLint)et);

	// TEXTURE->SAMPLERS (Per Step) \\

	// TEXTURE - DENSITY \\

	// Pass Density Grid - grid_data vector, data array to 3D Texture. 4Bytes (32Bit) Float Per Grid Density Value.
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_3D, tex_dens);
	// Use Linear (Trilinear Texture Filtering HC'd) 
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	GLfloat *ptr = (GLfloat*)f3obj->dens->grid_data->data();

	// 1 Channels (Red) Single Float for total Voxel/Texel per cell value -
	glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, (GLint) grid_size.x, (GLint) grid_size.y, (GLint) grid_size.z, 0, GL_RED, GL_FLOAT, ptr);

	// TEXTURE - VELOCITY \\

	// Pack X-Y-Z Velocity into R-G-B Components of single Velocity 3D Texture. 
	// Check Vel Grid Size (via grid_data flat array size) Actually == Passed GridSize to RenderObj. 
	assert(f3obj->vel->grid_data->size() == (grid_size.x * grid_size.y * grid_size.z)); // Passed Grid Size, should incl Edge Cells per Dim. 

	// Memory Layout - 
	// Per Voxel (3D Texel) - [Rx|Gy|Bz]. VelocityGrid Size (1D) * 3. Thus per voxel byte stride of 3 * sizeof (float)  

	std::size_t vel3dsize = f3obj->vel->grid_data->size() * 3;
	GLfloat *vel3D = new float[vel3dsize];

	// For Each Voxel, Write XYZ Vel Comps to RGB T/Voxel Comps - 
	for (std::size_t i = 0; i < f3obj->vel->grid_data->size(); i++)
	{
		// Get Velocity Components Of Current Cell - 
		float temp_U = f3obj->vel->getdata_x(i);
		float temp_V = f3obj->vel->getdata_y(i);
		float temp_W = f3obj->vel->getdata_z(i);

		// Per Voxel Components Set. 
		// XYZ -> RGB Components 
		for (std::size_t j = i; j < i+2; j++)
		{
			if (j % 3 == 0) vel3D[j] = temp_U; // R
			if (j % 3 == 1) vel3D[j] = temp_V; // G
			if (j % 3 == 2) vel3D[j] = temp_W; // B
		}

	}
	// Now Pass Single Flat 1D Array to 3D Velocity Texture - 

	// Velocity Texture - 
	glBindTexture(GL_TEXTURE_3D, 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_3D, tex_vel);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, grid_size.x, grid_size.y, grid_size.z, GL_RGB, GL_FLOAT, vel3D);

	// Delete Temp Arrays - 
	delete vel3D; vel3D = nullptr; 
}


// RenderObject_2D_OGL Shader Loader Implementation -
void renderobject_3D_OGL::render_loop(rend_state rs)
{
	// IF DBG Then Use While Loop Here - (As Assume NOT CALLED FROM INSIDE SOLVER LOOP (per step))
	// DBG Mode assumes RenderObject Setup/Debugging outside of a FluidSolver Instance eg from main for dbg sake. 
	if (rs == rend_state::RENDER_DEBUG)
	{
		while (!glfwWindowShouldClose(window_ptr))
		{
			// Render Loop 

			//glPolygonMode(GL_FRONT, GL_LINE);
			glClearColor(0.0f, 0.0f, 1.0f, 0.0f);

			// Active and Bind Textures. 
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, tex_dens);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, tex_vel);

			glUseProgram(shader_prog);

			// Draw Triangle in Render Loop
			glBindVertexArray(VAO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

			glfwSwapBuffers(window_ptr);
			glfwPollEvents();
		}

	}
	else // Assume Called INSIDE SOLVE LOOP. Thus RENDER_ACTIVE dbg state. 
	{
		//So Just do In Loop Operations (Because there already called within a (solve) loop) -	
		//Input Polling done within Solver Solve_Step Loop for more freedom of inputs vars. 

		// Render Loop \\

		glClearColor(0.0f, 0.0f, 1.0f, 0.0f);

		// Active and Bind Textures. (Multiple Texture/Units) - 
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, tex_dens);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_3D, tex_vel);

		// Call Shader Program. (Is also called in SolveStep on RenderObj Instance). 
		glUseProgram(shader_prog);

		this->print_GL_error(); // Check for GL Errors.
		
		// Draw Quad via Indexed Drawing - 
		glBindVertexArray(VAO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		glfwSwapInterval(0); // No VSync !
		glfwSwapBuffers(window_ptr);
		glfwPollEvents();
	}

}


// Temp Implementation to allow public acess to call render (From FluidSolver_2 SolveStep).
// FIX THIS HACKYNESS ! Why not Just Friend RenderObj to FluidSolver_3 ? 
void renderobject_3D_OGL::call_ren(rend_state rs)
{
	render_loop(rs);
}

void renderobject_3D_OGL::print_GL_error()
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
	{
		std::cout << "ERROR::RENDER_OBJECT_OGL:: " << err << "\n";
	}
}