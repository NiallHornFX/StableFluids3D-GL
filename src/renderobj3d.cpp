// Implementation of RenderObject3D
#include "renderobject3d.h"

// Project Headers
#include "fluidobj3d.h"
#include "assets.h"

// Std Headers
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

extern short verbose;

// Render Object Creation Classes Implementation - 

// RenderObject_3D - ABC Implementations - 

// Render Object 3D ABC Constructor
renderobject_3D::renderobject_3D(const char *api_name, int v_maj, int v_min, const vec2<int> &win_size, const vec3<int> &grd_size)
	: API_Name(api_name), ver_major(v_maj), ver_minor(v_min), window_size(win_size), grid_size(grd_size) {}

// RenderObject_2D_OGL (OpenGL API) Implementations - 

// RenderObject_2D_OGL Constructor - 
renderobject_3D_OGL::renderobject_3D_OGL(const char *api_name, int v_maj, int v_min, const vec2<int> &ws, const vec3<int> &gs, GLFWwindow *winptr)
	: renderobject_3D(api_name, v_maj, v_min, ws, gs), window_ptr(winptr),  // Initalize RenObj ABC Members Via Its Own Constructor. 
	cube_model(), cube_view(), cube_persp(), cam_target(0.0f, 0.0f, 0.0f), zoom(40.0f) // Ident Initalize Cube Transformation Matrix and Vector Members
{
	std::cout << "DBG::RenderObject_3D Created For Render API: " << api_name << " " << v_maj << "." << v_min << "\n \n";
	// Call Setup MFuncs

	cur_shader = 0;
	int sl_0 = shader_loader("shaders/Cube_BakeShader.vert", "shaders/Cube_BakeShader.frag");
	cur_shader = 1;
	int sl_1 = shader_loader("shaders/Quad_RayMarchingShader.vert", "shaders/Quad_RayMarchShader.frag");

	int vr = vertex_setup(); 

	// Gen 3D Textures -
	glGenTextures(1, &tex_dens);
	glGenTextures(1, &tex_vel);
	// Do Not Bind On Initalization.
}

// RenderObject_2D_OGL Destructor -
renderobject_3D_OGL::~renderobject_3D_OGL()
{
	// Dealloc Cube 

	if (cube_vert_shader || cube_vert_shader != NULL) {
		glDeleteShader(cube_vert_shader);
		cube_vert_shader = NULL;
	}

	if (cube_frag_shader || cube_frag_shader != NULL) {
		glDeleteShader(cube_frag_shader);
		cube_frag_shader = NULL;
	}

	if (cube_shader_prog || cube_shader_prog != NULL) {
		glDeleteProgram(cube_shader_prog);
		cube_shader_prog = NULL;
	}

	delete cube_vertices;  cube_vertices = nullptr;
	delete cube_edge_vertices; cube_edge_vertices = nullptr;
	delete cube_vert_shader_code; cube_vert_shader_code = nullptr;
	delete cube_frag_shader_code; cube_frag_shader_code = nullptr;

	// Dealloc Quad 

	if (quad_vert_shader || quad_vert_shader != NULL) {
		glDeleteShader(quad_vert_shader);
		quad_vert_shader = NULL;
	}

	if (quad_frag_shader || quad_frag_shader != NULL) {
		glDeleteShader(quad_frag_shader);
		quad_frag_shader = NULL;
	}

	if (quad_shader_prog || quad_shader_prog != NULL) {
		glDeleteProgram(quad_shader_prog);
		quad_shader_prog = NULL; 
	}

	delete quad_vertices; quad_vertices = nullptr; 
	delete quad_indices; quad_indices = nullptr; 
	delete quad_vert_shader_code; quad_vert_shader_code = nullptr; 
	delete quad_frag_shader_code; quad_frag_shader_code = nullptr; 
}

/* RenderObject_2D_OGL Vertex Setup Implementation - Setup Quad and Cube Buffers and Arrays Here */ 
int renderobject_3D_OGL::vertex_setup()
{
	// Cube Verts in assets.h copy included ptr to member ptr. 
	cube_vertices = cvert_n;
	cube_edge_vertices = cvert_edge; 

	// Vertex Screen Quad Setup - 
	quad_vertices = new GLfloat[12]
	{
		1.0f, 1.0f, 0.0f,    // 0 TR Vert
		1.0f, -1.0f, 0.0f,   // 1  BR Vert
	    -1.0f, -1.0f, 0.0f,  // 2 BL Vert
		-1.0f, 1.0f, 0.0f    // 3 TL Vert
	};

	quad_indices = new GLuint[6]
	{
		0,1,3, // Tri 1
		1,2,3  // Tri 2
	};

	// Cube Setup \\

	glGenVertexArrays(1, &Cube_VAO); glGenBuffers(1, &Cube_VBO);
	glBindVertexArray(Cube_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, Cube_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (36 * 6), cube_vertices, GL_STATIC_DRAW);
	// (float) XYZ-UVW (6 * sizeof(float) stride) 
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
	glEnableVertexAttribArray(0); // Enable VAO Attrib 0 - Postion
	// Vertex UV Attrib - (1) Start at float*3 offset of Postion Attribs. 
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(sizeof(float) * 3));
	glEnableVertexAttribArray(1);
	// UnBind
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Cube Edge Setup \\

	glGenVertexArrays(1, &Cube_Edge_VAO); glGenBuffers(1, &Cube_Edge_VBO);
	glBindVertexArray(Cube_Edge_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, Cube_Edge_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (24 * 6), cube_edge_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
	glEnableVertexAttribArray(0); // Attrib 0 - Postion
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(sizeof(float) * 3));
	glEnableVertexAttribArray(1); // Attrib 1 - Colour
	// UnBind
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Full Screen Quad Setup \\ 

	glGenVertexArrays(1, &Quad_VAO);
	glGenBuffers(1, &Quad_VBO);
	glBindVertexArray(Quad_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, Quad_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12, quad_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0); // VertPos
	glEnableVertexAttribArray(0);
	// UnBind
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Quad Element Index Buffer Setup - 
	glGenBuffers(1, &Quad_EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Quad_EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 6, quad_indices, GL_STATIC_DRAW);
	// UnBind
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// QUAD Texture Sampler Units - 
	glUseProgram(quad_shader_prog);
	// Set Quad Shader 2D Texture Sampler Units (Cube Front Face and Back Face 0,1) -
	glUniform1i(glGetUniformLocation(quad_shader_prog, "cf_tex"), 0); 
	glUniform1i(glGetUniformLocation(quad_shader_prog, "cb_tex"), 1); 
	// Set Quad Shader 3D Texture Sampler Units (Density, Velocity Grid Textures 2,3) - 
	glUniform1i(glGetUniformLocation(quad_shader_prog, "dens_tex"), 2); 
	glUniform1i(glGetUniformLocation(quad_shader_prog, "vel_tex"), 3); 
	glUseProgram(0);

	// Inital RenderState 
	inital_renderstate();

	return 0;
}

// RenderObject_3D_OGL Shader Compilation Checker MFunc Implementation - 
void renderobject_3D_OGL::shader_checkCompile(const char *type, int shader)
{
	assert(type == (const char*) "vertex" || type == (const char*) "fragment");

	const int len = 512; // length of error log char array. 
	int sucess_v, sucess_f;
	char err_log_v[len], err_log_f[len];

	// Vertex Shader Check - 
	if (strcmp(type, "vertex") == 0)
	{
		GLuint cur_vert_shader;
		if (shader == 0) { cur_vert_shader = cube_vert_shader; } else if (shader == 1) { cur_vert_shader = quad_vert_shader; }

		glGetShaderiv(cur_vert_shader, GL_COMPILE_STATUS, &sucess_v);
		if (!sucess_v)
		{
			glGetShaderInfoLog(cur_vert_shader, len, NULL, err_log_v);
			std::cerr << "ERR::VERTEX SHADER " << shader << " - " << cur_vert_shader << " COMPILE FAILED \n " << err_log_v << std::endl; 
			std::terminate();
		}
		else if (verbose == 1)
		{
			std::cout << "DBG::VERTEX SHADER COMPILE SUCEEDED \n";
		}
	}
	
	// Fragment Shader Check - 
	if (strcmp(type, "fragment") == 0)
	{
		GLuint cur_frag_shader;
		if (shader == 0) { cur_frag_shader = cube_frag_shader; } else if (shader == 1) { cur_frag_shader = quad_frag_shader; }

		glGetShaderiv(cur_frag_shader, GL_COMPILE_STATUS, &sucess_f);
		if (!sucess_f)
		{
			glGetShaderInfoLog(cur_frag_shader, len, NULL, err_log_f);
			std::cerr << "ERR::FRAGMENT SHADER " << shader << " - " << cur_frag_shader << " COMPILE FAILED \n " << err_log_f << std::endl;
			std::terminate();
		}
		else if (verbose == 1)
		{
			std::cout << "DBG::FRAGMENT SHADER COMPILE SUCEEDED \n";
		}
	}
}

// RenderObject_3D_OGL ShaderProgram Linker Checker MFunc Implementation - 
void renderobject_3D_OGL::shader_checkLink(int shader)
{
	GLuint cur_shader_prog;
	if (shader == 0) { cur_shader_prog = cube_shader_prog; } else if (shader == 1) { cur_shader_prog = quad_shader_prog; }
	std::cout << "SHADER " << shader << " PROG " << cur_shader_prog << "\n";
	int sucess;
	const int len = 512;
	char err_log[len];

	glGetProgramiv(cur_shader_prog, GL_LINK_STATUS, &sucess);

	if (!sucess)
	{
		glGetProgramInfoLog(cur_shader_prog, len, NULL, err_log);
		std::cerr << "ERR:SHADER-PROGRAM: " << shader << " - " << cur_shader_prog << " LINKAGE_FAILED" << std::endl;
		std::cerr << err_log << std::endl;
		std::terminate();
	}
	else if (verbose == 1)
	{
		std::cout << "DBG::SHADER:PROGRAM " << cur_shader_prog << " LINKAGE SUCEEDED" << std::endl;
	}
}

// RenderObject_3D_OGL Shader Loader Implementation - 
int renderobject_3D_OGL::shader_loader(const char *vert_path, const char *frag_path)
{
	std::ifstream vert_shader_load, frag_shader_load; 
	std::stringstream v_shad_buf, f_shad_buf; 
	std::string temp_v, temp_f; 

	vert_shader_load.exceptions(std::ios::badbit | std::ios::failbit);
	frag_shader_load.exceptions(std::ios::badbit | std::ios::failbit);

	/* Shader Load And Write to Code Buffer. (Probs could of done this with strcpy directly much easier !) 
	 No Parsing, (For now) just check compile and linkage on shader creation below. */ 
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
		std::terminate();
	}

	if (cur_shader == 0) // Cube 
	{
		// Fill Shader Code Member CStrs. 
		cube_vert_shader_code = temp_v.c_str();
		cube_frag_shader_code = temp_f.c_str();

		// Cube Shaders Compilation & Link - 
		cube_vert_shader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(cube_vert_shader, 1, &cube_vert_shader_code, NULL);
		glCompileShader(cube_vert_shader);

		cube_frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(cube_frag_shader, 1, &cube_frag_shader_code, NULL);
		glCompileShader(cube_frag_shader);

		shader_checkCompile("vertex", 0); shader_checkCompile("fragment", 0);

		// Cube Shader Program - 
		cube_shader_prog = glCreateProgram();
		glAttachShader(cube_shader_prog, cube_vert_shader); glAttachShader(cube_shader_prog, cube_frag_shader);
		glLinkProgram(cube_shader_prog);

		shader_checkLink(0);
	} 
	else if (cur_shader == 1) // Quad
	{
		// Fill Shader Code Member CStrs. 
		quad_vert_shader_code = temp_v.c_str();
		quad_frag_shader_code = temp_f.c_str();

		// Quad Shaders Compilation & Link - 
		quad_vert_shader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(quad_vert_shader, 1, &quad_vert_shader_code, NULL);
		glCompileShader(quad_vert_shader);

		quad_frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(quad_frag_shader, 1, &quad_frag_shader_code, NULL);
		glCompileShader(quad_frag_shader);

		shader_checkCompile("vertex", 1); shader_checkCompile("fragment", 1);

		// Quad Shader Program - 
		quad_shader_prog = glCreateProgram();
		glAttachShader(quad_shader_prog, quad_vert_shader); glAttachShader(quad_shader_prog, quad_frag_shader);
		glLinkProgram(quad_shader_prog);

		shader_checkLink(1);
	}

	return 0; 
}

// Call Per Frame to Convert and Update Grids to 3D Textures and utility uniforms for Rendering.
void renderobject_3D_OGL::shader_pipe(fluidobj_3d *f3obj)
{
	// Shader Pipe is Passing to Quad Shader Only.
	glUseProgram(quad_shader_prog);

	// TEXTURE->SAMPLERS (Per Step) \\
	// TEXTURE - DENSITY \\

	glBindTexture(GL_TEXTURE_3D, 0);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_3D, tex_dens);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	GLfloat *ptr = (GLfloat*)f3obj->dens->grid_data->data();

	// 1 Channels (Red) Single Float for total Voxel/Texel per cell value -
	glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, (GLint) grid_size.x, (GLint) grid_size.y, (GLint) grid_size.z, 0, GL_RED, GL_FLOAT, ptr);

	// TEXTURE - VELOCITY \\

	// Pack X-Y-Z Velocity into R-G-B Components of single Velocity 3D Texture.  
	assert(f3obj->vel->grid_data->size() == (grid_size.x * grid_size.y * grid_size.z)); 

	// Per Voxel (3D Texel) - [Rx|Gy|Bz]. VelocityGrid Size (1D) * 3. Thus per voxel byte stride of 3 * sizeof (float)  

	std::size_t vel3dsize = f3obj->vel->grid_data->size() * 3;
	auto vel3D = std::make_unique<GLfloat[]>(vel3dsize); 

	// Merge XYZ Vel Components to 3D Scalar RGB Componets Per Cell. 
	for (int i = 0, ch_idx = 0; i < f3obj->t_s; i++, ch_idx += 3)
	{
		vel3D[ch_idx] = f3obj->vel->getdata_x(i); // v.x -> R
		vel3D[ch_idx + 1] = f3obj->vel->getdata_y(i); // v.y -> G
		vel3D[ch_idx + 2] = f3obj->vel->getdata_z(i); // v.z -> B
	}

	// Flat 1D Array to 3D Velocity Texture -  
	glBindTexture(GL_TEXTURE_3D, 0);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_3D, tex_vel);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, (GLint)grid_size.x, (GLint)grid_size.y, (GLint)grid_size.z, 0, GL_RGB, GL_FLOAT, vel3D.get()); 

	glBindTexture(GL_TEXTURE_3D, 0);
	glUseProgram(0);

}

// Need a better way to pass this, most likely through shader_pipe().
void renderobject_3D_OGL::get_input(const vec2<float> &m)
{
	glUseProgram(quad_shader_prog);
	glUniform1f(glGetUniformLocation(quad_shader_prog, "mx"), m.x);
	glUniform1f(glGetUniformLocation(quad_shader_prog, "my"), m.y);
	glUseProgram(0);
	std::cout << std::fixed << "MOUSE POS = " << m.x << "  " << m.y << "\n";
}

// RenderObject_3D_OGL Cube Setup - Setup Cube Transforms (Initalize, and then cube_update() in RLoop?)
void renderobject_3D_OGL::cube_setup()
{
	// Inital Cube Transform Setup to pass to GPU \\ 

	// Model-World Matrix - 
	cube_model.rotate(vec3<float>(0.0f, 1.0f, 0.0f), -0.1f);
	//cube_model.rotate(vec3<float>(1.0f, 1.0f, 0.0f), 0.1f);
	cube_model.label = "Cube Model";
	cube_model.print_mat();

	// View Matrix - 
	cube_view.translate(vec3<float>(0.0f, 0.0f, -3.0f)); // No LA Yet. Just Move Back on -z (ie cam "moved" along +z)
	cube_view.label = "Cube View";
	cube_view.print_mat();

	// Perspective Proj Matrix - 
	cube_persp = matrix_4x4<float>::make_perspective(40.0f, 1.0f, 0.1f, 100.0f);
	cube_persp.label = "Cube Persp";
	cube_persp.print_mat();

	glUseProgram(cube_shader_prog);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "model"), 1, GL_FALSE, cube_model.comp);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "view"), 1, GL_TRUE, cube_view.comp);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "persp"), 1, GL_FALSE, cube_persp.comp);
	glUseProgram(0);
}

// Update Transformation Matrices in Render loop Testing, Need to add Parms. 
void renderobject_3D_OGL::cube_update()
{
	//cube_model.rotate(vec3<float>(0.0f, 1.0f, 0.0f), matrix_4x4<float>::degtoRad(-(2.0 * dt)));
	float ang = matrix_4x4<float>::degtoRad( (std::sinf(t1 * 1.5f) * -15.0f) * dt); // Some Fake Camera/Grid Motion.
	ang = t1 * 0.05f; //* (1.0f * dt); 
	std::cout << "DEBUG::Rotation Angle = " << t1 << "\n"; 
	//cube_model.rotate(vec3<float>(0.0f, 1.0f, 0.0f), -ang);
	cube_model = matrix_4x4<float>::make_rotate(vec3<float>(0.0f, 1.0f, 0.0f), -ang); // New each frame, oppose to product of prev.

	// Testing Camera Translation like beahviour. 
	if (glfwGetKey(window_ptr, GLFW_KEY_W) == GLFW_PRESS)
	{
		cube_view.translate(vec3<float>(0.0f, 0.0f, 1.0f * dt));
	}
	if (glfwGetKey(window_ptr, GLFW_KEY_S) == GLFW_PRESS)
	{
		cube_view.translate(vec3<float>(0.0f, 0.0f, -1.0f * dt));
	}

	// If FOV Changed Update Persp Mat. 
	if (glfwGetKey(window_ptr, GLFW_KEY_UP) == GLFW_PRESS)
	{
		zoom -= 0.25f; 
	}
	if (glfwGetKey(window_ptr, GLFW_KEY_DOWN) == GLFW_PRESS)
	{
		zoom += 0.25f;
	}
	cube_persp = matrix_4x4<float>::make_perspective(zoom, 1.0f, 0.1f, 100.0f);

	// Update Uniforms - 
	glUseProgram(cube_shader_prog);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "model"), 1, GL_FALSE, cube_model.comp);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "view"), 1, GL_TRUE, cube_view.comp);
	glUniformMatrix4fv(glGetUniformLocation(cube_shader_prog, "persp"), 1, GL_FALSE, cube_persp.comp);
	glUseProgram(0);
}

// RenderObject_3D_OGL Cube FrameBuffer Setup - For Baking Rasterized Cube to - 
void renderobject_3D_OGL::cube_fbo_setup()
{
	// Gen FBOs - MultiSampled and Standard.  
	glGenFramebuffers(1, &Cube_FBO_s); glGenFramebuffers(1, &Cube_FBO_ms);

	// Gen RBO - MultiSampled.
	glGenRenderbuffers(1, &Cube_RBO);
	glBindRenderbuffer(GL_RENDERBUFFER, Cube_RBO);
	glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH24_STENCIL8, window_size.x, window_size.y);
	//glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, window_size.x, window_size.y); nMS
	glBindRenderbuffer(GL_RENDERBUFFER, 0); // Unbind RBO State. 

	// Gen MultiSampled Texture
	glGenTextures(1, &tex_ms);
	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, tex_ms);
	glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, GL_RGB, window_size.x, window_size.y, GL_TRUE);
	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

	// Generate Standard Textures (Front and Back Final Output) -
	// Front 
	glGenTextures(1, &tex_CFront_s);
	glBindTexture(GL_TEXTURE_2D, tex_CFront_s);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, window_size.x, window_size.y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);
	// Back
	glGenTextures(1, &tex_CBack_s);
	glBindTexture(GL_TEXTURE_2D, tex_CBack_s);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, window_size.x, window_size.y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	// Bind MS FBO and Attach MSTexture and RBO (DepthStencill nonMS) to it.  
	glBindFramebuffer(GL_FRAMEBUFFER, Cube_FBO_ms);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, tex_ms, 0); // Tex Colour (Will be Switched in RLoop)
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, Cube_RBO); // RBO DepthStenicll. NNTB
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) std::cerr << "ERR::FrameBuffer Not Complete. Check Attachments. \n"; 

	// Unbind remaing bound state.
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

}

// Bind and Clear Set FBO Ready for drawing to. Note Bind/Clear of FBO_S is done within fbo_attachAndblit before each MS-S Blit. 
void renderobject_3D_OGL::bindclear_fbo(use_fbo mode)
{
	if (mode == FBO_CUBE_MS)
	{
		// Cube FBO - With Depth Test. 
		glBindFramebuffer(GL_FRAMEBUFFER, Cube_FBO_ms);
		glEnable(GL_DEPTH_TEST);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}
	else if (mode == FBO_DEFAULT)
	{
		// Default FBO - No Depth Test (Quad Draw). 
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glDisable(GL_DEPTH_TEST);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
	}

}

/* RenderObject_3D_OGL - CubeFBO_AttachandBlit
	Attach - Switch FBO_Standard Cube Texture Attachments to Blit to. Colour Only. No Standard RBO. 
	Blit - Clear FBO_Standard Colour and Blit from FBO_MultiSample to it. */
void renderobject_3D_OGL::cube_fbo_attachAndblit(use_cube tex)
{
	if (tex == CUBE_FRONT) // FBO --> Cube Front Texture Attachment
	{
		glBindFramebuffer(GL_FRAMEBUFFER, Cube_FBO_s);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_CFront_s, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
	else if (tex == CUBE_BACK) // FBO --> Cube Back Texture Attachemnt
	{
		glBindFramebuffer(GL_FRAMEBUFFER, Cube_FBO_s);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_CBack_s, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	// Clear FBO_S Before Blitting - 
	glBindFramebuffer(GL_FRAMEBUFFER, Cube_FBO_s);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// Blit to FBO_S and its current colour texture attachment - 
	glBindFramebuffer(GL_READ_FRAMEBUFFER, Cube_FBO_ms);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, Cube_FBO_s);
	glBlitFramebuffer(0, 0, window_size.x, window_size.y, 0, 0, window_size.x, window_size.y, GL_COLOR_BUFFER_BIT, GL_LINEAR);
	// Clear all bound FBO State. 
	glBindFramebuffer(GL_READ_FRAMEBUFFER, 0); glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// Only Call this Once on RenderObject Initalization, After Vertex/Shader setup in RenderObject Ctor.
void renderobject_3D_OGL::inital_renderstate()
{
	cube_setup(); // Intial Cube Transform Setup
	cube_fbo_setup(); // Inital FBO Setup 
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glEnable(GL_MULTISAMPLE);
}

void renderobject_3D_OGL::get_FPS()
{
	step++;
	t1 = glfwGetTime();
	dt = t1 - t0;
	t0 = t1;
	if (step % 1 == 0) std::cout << std::fixed << "DEBUG::RENDER_OBJ FPS = " << 1.0f / dt << " FPS" << "\n"; // !TOD Logger Class / Utilsh. 
}

// RenderObject_3D_OGL Shader Loader Implementation - Debug State For Out Of Solver Testing Only. 
void renderobject_3D_OGL::render_loop(rend_state rs)
{
	if (rs == rend_state::RENDER_DEBUG)
	{
		while (!glfwWindowShouldClose(window_ptr))
		{
			// PRE OP 
			get_FPS();

			// CUBE \\ 
			cube_update(); // Transforms
			glUseProgram(cube_shader_prog);

			// CUBE FRONT 
			bindclear_fbo(FBO_CUBE_MS);
			// Draw Cube Front -
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK); 
			glBindVertexArray(Cube_VAO);
			glDrawArrays(GL_TRIANGLES, 0, 36);
			// Blit to Standard FrameBuffer Texture (Front) - 
			cube_fbo_attachAndblit(CUBE_FRONT);

			// CUBE BACK 
			bindclear_fbo(FBO_CUBE_MS); 
			// Draw Cube Back -
			glEnable(GL_CULL_FACE);
			glCullFace(GL_FRONT);
			glBindVertexArray(Cube_VAO);
			glDrawArrays(GL_TRIANGLES, 0, 36);
			// Blit to Standard FrameBuffer Texture (Back) - 
			cube_fbo_attachAndblit(CUBE_BACK);

			// SCREEN QUAD 
			bindclear_fbo(FBO_DEFAULT);
			// Reset Render State
			glDisable(GL_CULL_FACE);
			glUseProgram(0);
			glBindVertexArray(0);

			// Render Screen Quad - 
			glUseProgram(quad_shader_prog);

			// Active and Bind Textures. 
			// 2D - 
			glActiveTexture(GL_TEXTURE0); // Cube Front. 
			glBindTexture(GL_TEXTURE_2D, tex_CFront_s);
			glActiveTexture(GL_TEXTURE1); // Cube Back. 
			glBindTexture(GL_TEXTURE_2D, tex_CBack_s);
			// 3D - 
			glActiveTexture(GL_TEXTURE2); // Density Grid
			glBindTexture(GL_TEXTURE_3D, tex_dens);
			glActiveTexture(GL_TEXTURE3); // Velocity Grid
			glBindTexture(GL_TEXTURE_3D, tex_vel);

			// Draw Full Screen Quad to Default FrameBuffer. 
			glBindVertexArray(Quad_VAO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Quad_EBO);
			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

			// POST OP \\ 

			// Clear Render State. 
			glUseProgram(0);
			glBindVertexArray(0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

			// Swap and Poll - 
			glfwSwapInterval(0); // Disable vsync. 
			glfwSwapBuffers(window_ptr);
			glfwPollEvents();
		}

	}
	else if (rs == rend_state::RENDER_ACTIVE)
	{
		// PRE OP 
		get_FPS();

		// CUBE \\ 
		cube_update(); // Transforms
		glUseProgram(cube_shader_prog);

		// CUBE FRONT 
		bindclear_fbo(FBO_CUBE_MS);
		// Draw Cube Front -
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glBindVertexArray(Cube_VAO);
		glDrawArrays(GL_TRIANGLES, 0, 36);
		// Blit to Standard FrameBuffer Texture (Front) - 
		cube_fbo_attachAndblit(CUBE_FRONT);

		// CUBE BACK 
		bindclear_fbo(FBO_CUBE_MS);
		// Draw Cube Back -
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glBindVertexArray(Cube_VAO);
		glDrawArrays(GL_TRIANGLES, 0, 36);
		// Blit to Standard FrameBuffer Texture (Back) - 
		cube_fbo_attachAndblit(CUBE_BACK);

		// SCREEN QUAD 
		bindclear_fbo(FBO_DEFAULT);
		// Reset Render State
		glDisable(GL_CULL_FACE);
		glUseProgram(0);
		glBindVertexArray(0);

		// Render Screen Quad - 
		glUseProgram(quad_shader_prog);

		// Active and Bind Textures. 
		// 2D - 
		glActiveTexture(GL_TEXTURE0); // Cube Front. 
		glBindTexture(GL_TEXTURE_2D, tex_CFront_s);
		glActiveTexture(GL_TEXTURE1); // Cube Back. 
		glBindTexture(GL_TEXTURE_2D, tex_CBack_s);
		// 3D - 
		glActiveTexture(GL_TEXTURE2); // Density Grid
		glBindTexture(GL_TEXTURE_3D, tex_dens);
		glActiveTexture(GL_TEXTURE3); // Velocity Grid
		glBindTexture(GL_TEXTURE_3D, tex_vel);

		// Draw Full Screen Quad to Default FrameBuffer. 
		glBindVertexArray(Quad_VAO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Quad_EBO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		// Draw Cube Edges to Default FrameBuffer atop. 
		glUseProgram(0);
		glUseProgram(cube_shader_prog);
		glEnable(GL_MULTISAMPLE);
		glBindVertexArray(Cube_Edge_VAO);
		glDrawArrays(GL_LINES, 0, 24);

		// POST OP \\ 

		// Clear Render State. 
		glUseProgram(0);
		glBindVertexArray(0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		// Swap and Poll - 
		glfwSwapInterval(0); // Disable vsync. 
		glfwSwapBuffers(window_ptr);
		glfwPollEvents();
	}

}

void renderobject_3D_OGL::print_GL_error()
{
	GLenum err;
	int i = 0;
	while ((err = glGetError()) != GL_NO_ERROR && i < 5)
	{
		std::cout << "ERROR::RENDER_OBJECT_OGL:: " << err << "\n";
		i++; 
	}
}