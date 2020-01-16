//#include "renderobject.h"
//#include "fluidobj.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

extern short verbose;

// Render Object Creation Classes Implementation - 

// RenderObject_2D - ABC Implementations - 

// Render Object 2D ABC Constructor
renderobject_2D::renderobject_2D(const char *api_name, int v_maj, int v_min, int win_x, int win_y, int rmode)
	: API_Name(api_name), ver_major(v_maj), ver_minor(v_min), winsize_x(win_x), winsize_y(win_y), rendermode(rmode) {}


// RenderObject_2D_OGL (OpenGL API) Implementations - 

// RenderObject_2D_OGL Constructor - 
renderobject_2D_OGL::renderobject_2D_OGL(const char *api_name, int v_maj, int v_min, int wx, int wy, GLFWwindow *winptr, int rmode)
	: renderobject_2D(api_name, v_maj, v_min, wx, wy, rmode), window_ptr(winptr)  // Initalize ABC Members Via Its Own Constructor. 
{
	std::cout << "DBG::RenderObject_2D Created For Render API: " << api_name << " " << v_maj << "." << v_min << "\n \n";
	// Call Setup MFuncs

	int v_ret = vertex_setup(); 
	int s_ret = shader_loader("render2d_ogl_vertShader.vert", "render2d_ogl_fragShader.frag" ); // Hard Coded Shader Paths For now. 

	glGenTextures(1, &tex_dens);
	glGenTextures(1, &tex_vel_u);
	glGenTextures(1, &tex_vel_v);
	glGenTextures(1, &tex_c);
	glGenTextures(1, &tex_vc_u);
	glGenTextures(1, &tex_vc_v);
	glGenTextures(1, &tex_img_rgb);
	glGenTextures(1, &tex_preprojvel_u);
	glGenTextures(1, &tex_preprojvel_v);
	// Dont Bind Any yet.
}

// RenderObject_2D_OGL Destructor -
renderobject_2D_OGL::~renderobject_2D_OGL()
{
	if (vert_shader != NULL) {
		glDeleteShader(vert_shader);
		vert_shader = NULL;
	}

	if (frag_shader != NULL) {
		glDeleteShader(frag_shader);
		frag_shader = NULL;
	}

	if (shader_prog != NULL) {
		glDeleteProgram(shader_prog);
		shader_prog = NULL; 
	}

	delete vertices; vertices = nullptr; 
	delete indices; indices = nullptr; 

	delete vert_shader_code; vert_shader_code = nullptr; 
	delete frag_shader_code; frag_shader_code = nullptr; 
}

// RenderObject_2D_OGL Vertex Setup Implementation, For Quad and Indexed Vertex Drawing - 
int renderobject_2D_OGL::vertex_setup()
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

	// Vertex Array/Buffer Setup - 
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

// RenderObject_2D_OGL Shader Compilation Checker MFunc Implementation - 
void renderobject_2D_OGL::shader_checkCompile(const char *type)
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

// RenderObject_2D_OGL ShaderProgram Linker Checker MFunc Implementation - 
void renderobject_2D_OGL::shader_checkLink()
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

// RenderObject_2D_OGL Shader Loader Implementation - 
int renderobject_2D_OGL::shader_loader(const char *vert_path, const char *frag_path)
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

	// Dealloc Shaders Here vs Destruc? 
	//glDeleteShader(vert_shader); vert_shader = NULL; 
	//glDeleteShader(frag_shader); frag_shader = NULL;

	return 0; 
}

void renderobject_2D_OGL::shader_pipe(fluidobj_2d *f2obj)
{
	// UNIFORM CONSTANTS (Per Step) \\

	// Pass Window Size (Single_Dimension / N + EdgeSize) to GL Uniform. // Could be done after Shader Compile (Once.?)
	glUniform1i(glGetUniformLocation(shader_prog, "N_Size"), winsize_x); // Assuming Square Grid/Window so just pass X size. 
	//std::cout << "DEBUG ***********    " << winsize_x << "\n \n";

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

	// Pass Density Grid - grid_data vector, data array to 2D Texture. 4Bytes (32Bit) Float Per Grid Density Value.
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_dens);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	GLfloat *ptr = (GLfloat*)f2obj->dens->grid_data->data(); // Dont need this cast to GLFloat, but good for explicit. 
	// 1 Channel (Red) Single Float for total texel per cell value -
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, ptr);

	// TEXTURE - VELOCITY \\

	// Pass Vel Grid Per Component - grid_data vector, data array to 2D Texture. 4Bytes (32Bit) Float Per Grid Velocity u/v value.
	// Need to Sepreate Vel vec2 into u and v Float Component Arrays. 
	// Pass Vel Components to Temp Per Component Array to Pass to Vel_U Texture (This is not ideal).  Ideally Split within gridobject method.
	GLfloat *temp_u = new float[int(f2obj->vel->grid_data->size())];
	GLfloat *temp_v = new float[int(f2obj->vel->grid_data->size())];
	// Parrelize ? 
	for (int i = 0; i < f2obj->vel->grid_data->size(); i++)
	{
		temp_u[i] = std::fabsf(f2obj->vel->grid_data->at(i).x);
		temp_v[i] = std::fabsf(f2obj->vel->grid_data->at(i).y);
	}
	// Velocity X U Component Texture -
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, tex_vel_u);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_u);

	// Velocity Y V Component Texture -
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, tex_vel_v);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_v);

	// Delete Temp Vel Component Grids 
	delete[] temp_u; temp_u = nullptr;
	delete[] temp_v; temp_v = nullptr; 

	// TEXTURE - COLLISION \\

	// Render Col (Collider) Grid (if passed) - 
	if (f2obj->col)
	{
		glBindTexture(GL_TEXTURE_2D, 0);
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, tex_c);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		GLfloat *col_cast = (GLfloat*)f2obj->col->grid_data->data();
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, col_cast);
	}

	// TEXTURE - VORTICITY \\

	// Render Vorticity Confinement Grid (Only if Exsits) - 
	// Like Velocity Split u,v VortConfine Components into 32bit Float Single Channel (GL_Red) Textures to pass to GPU via GLTexImage2D.
	
	if (f2obj->vc)
	{
		// Split Vort Confine vec2 Grid Components into Seperate u,v float arrays, for each u,v tex. 
		GLfloat *temp_vc_u = new float[int(f2obj->vc->grid_data->size())];
		GLfloat *temp_vc_v = new float[int(f2obj->vc->grid_data->size())];
		for (int i = 0; i < f2obj->vc->grid_data->size(); i++)
		{
			temp_vc_u[i] = std::fabsf(f2obj->vc->grid_data->at(i).x);
			temp_vc_v[i] = std::fabsf(f2obj->vc->grid_data->at(i).y);
		}

		// Vorticty Confinement X U Component Texture - 
		glBindTexture(GL_TEXTURE_2D, 0);
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, tex_vc_u);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_vc_u);

		// Vorticty Confinement Y V Component Texture - 
		glBindTexture(GL_TEXTURE_2D, 0);
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, tex_vc_v);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_vc_v);

		// Delete Temp Vel Component Grids 
		delete[] temp_vc_u;
		delete[] temp_vc_v;
	}

	// TEXTURE - RGB \\

	// RGB Image (Colour) Grid Texture - 
	glActiveTexture(GL_TEXTURE6);
	glBindTexture(GL_TEXTURE_2D, tex_img_rgb);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	// Merge R,G,B Grids into to One Float Array (GL_RGB). 
	const unsigned int merge_size = f2obj->c_R->grid_data->size() + f2obj->c_G->grid_data->size() + f2obj->c_B->grid_data->size(); 
	GLfloat *RGBMerge = new GLfloat[merge_size];
	assert(sizeof(GLfloat) == sizeof(float)); // Incase diff precision, will have bounds issues (grid vs tex stride).

	// Assume Colour Grids are same dim (they should be) - 
	assert(f2obj->c_R->grid_data->size() == f2obj->c_G->grid_data->size());
	const unsigned int grid_size = f2obj->c_R->grid_data->size(); 

	// For All Cells Get RGB Values, and write to Merged {R-G-B} (*3) Per Pixel Array to pass to Texture. 
	// Do All Cells, assume Edge Cells do not have colour from sim itself, oppose to excluding them here. 
	int ch_idx = 0; 
	for (int i = 0; i < grid_size; i++)
	{
		// Get Colour Of Pixel From Colour Cur Grid Cell
		GLfloat R_v = (GLfloat)f2obj->c_R->getdata(i) / 255.0f;
		GLfloat G_v = (GLfloat)f2obj->c_G->getdata(i) / 255.0f;
		GLfloat B_v = (GLfloat)f2obj->c_B->getdata(i) / 255.0f;

		// Write to Mered Texture RGB Array - 
		RGBMerge[ch_idx] = R_v; 
		RGBMerge[ch_idx+1] = G_v;
		RGBMerge[ch_idx+2] = B_v;

		ch_idx += 3; // Stride Offset Num of Channels Per Pixel. 
	}

	// Pass Data To GPU Texture -
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, winsize_x, winsize_y, 0, GL_RGB, GL_FLOAT, RGBMerge);

	// Clear CPU Resources 
	delete[] RGBMerge; RGBMerge = nullptr; 

	// TEXTURE - DEBUG - PRE-PROJECTED VELOCITY \\ 

	// Draw PreProjected Velocity for debug purposes, and visulization of residual of pressure solve and HH Decomp. 
	// Split PreProjVel vec2 components into Temp Arrays to pass to PPVel U,V Textures -
	GLfloat *temp_pp_u = new float[int(f2obj->preproj_vel->grid_data->size())];
	GLfloat *temp_pp_v = new float[int(f2obj->preproj_vel->grid_data->size())];

	for (int i = 0; i < f2obj->vel->grid_data->size(); i++)
	{
		temp_pp_u[i] = std::fabsf(f2obj->preproj_vel->grid_data->at(i).x);
		temp_pp_v[i] = std::fabsf(f2obj->preproj_vel->grid_data->at(i).y);
	}
	// Velocity X U Component Texture -
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D, tex_preprojvel_u);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_pp_u);

	// Velocity Y V Component Texture -
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE8);
	glBindTexture(GL_TEXTURE_2D, tex_preprojvel_v);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, winsize_x, winsize_y, 0, GL_RED, GL_FLOAT, temp_pp_v);

	// Delete Temp Vel Component Grids 
	delete[] temp_pp_u; temp_pp_u = nullptr;
	delete[] temp_pp_v; temp_pp_v = nullptr;

}


// RenderObject_2D_OGL Shader Loader Implementation -
void renderobject_2D_OGL::render_loop(rend_state dbg)
{
	// IF DBG Then Use While Loop Here - (As Assume NOT CALLED FROM INSIDE SOLVER LOOP (per step))
	// DBG Mode assumes RenderObject Setup/Debugging outside of a FluidSolver Instance eg from main for dbg sake. 
	// Now Using Enum to identify these states, oppose to bool true/false which is not very clear. 
	if (dbg == rend_state::RENDER_DEBUG)
	{
		// Texture Debug Code, Ignore...
		int count = 0;
		glBindTexture(GL_TEXTURE_2D, tex_dens);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		unsigned char *barr = new unsigned char[66564 * 3];
		//memset(barr, 1000.0f, 66564 * 3);
		for (int i = 0; i < int(66564 * 3); i++)
		{
			barr[i] = float(abs(sin(i))) * 250.0f;
		}
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 258, 258, 0, GL_RGB, GL_UNSIGNED_BYTE, barr);
		glUniform1i(glGetUniformLocation(shader_prog, "d_tex"), 0); // Enfroce Texture Unit 0. 
																	//delete barr;
		while (!glfwWindowShouldClose(window_ptr))
		{
			// Do Input Polling in hear for now? 
			// Render Loop 
			std::cout << count << "\n";
			count++;

			//glPolygonMode(GL_FRONT, GL_LINE);
			glClearColor(0.0f, 0.0f, 1.0f, 0.0f);

			// Active and Bind Textures. 
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, tex_dens);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, tex_vel_u);

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
		// Do Input Polling in hear for now? No. Input Polling done within Solver Solve_Step Loop for more freedom of inputs vars. 

		//if (glfwWindowShouldClose(window_ptr)) return; // Kill RenderStep if GLFW WindowClose. (Wont Kill Solver).

		// Render Loop \\

		//glPolygonMode(GL_FRONT, GL_LINE);
		glClearColor(0.0f, 0.0f, 1.0f, 0.0f);

		// Active and Bind Textures. (Multiple Texture/Units) - 
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, tex_dens);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, tex_vel_u);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, tex_vel_v);

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, tex_c);

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, tex_vc_u);
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, tex_vc_v);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, tex_img_rgb);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, tex_preprojvel_u);
		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, tex_preprojvel_v);

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
// FIX THIS HACKYNESS ! Why not Just Friend RenderObj to FluidSolver_2 ? 
void renderobject_2D_OGL::call_ren(rend_state dbg)
{
	render_loop(dbg);
}

void renderobject_2D_OGL::print_GL_error()
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
	{
		std::cout << "ERROR::RENDER_OBJECT_OGL:: " << err << "\n";
	}
}