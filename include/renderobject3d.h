#ifndef RENDEROBJECT_H
#define RENDEROBJECT_H

// Prevent Glew and GLFW multiple inclusion/implmentation defintion. 
#ifndef GLEW_STATIC
#define GLEW_STATIC  
#include <GLEW\glew.h>
#include <GLFW\glfw3.h>
#endif 

#include "grids3d.h"

// Forward Decls - 
class fluidobj_3d; 
class fluidsolver_3;

// Render State Enum 
enum class rend_state
{
	RENDER_DEBUG,
	RENDER_ACTIVE
};

enum tex_interp
{
	NEAREST = 0,
	TRILINEAR
};


// Abstract Base Class Interface for API Specfic Render Object. 
class renderobject_3D
{
// FCs
friend class fluidsolver_3;
public:
	renderobject_3D(const char *api_name, int v_maj, int v_min, const vec2<int> &win_size, const vec3<int> &grd_size);
	virtual ~renderobject_3D() = default;

	virtual void shader_pipe(fluidobj_3d *f3obj) = 0;
	virtual void render_loop(rend_state rs) = 0;

protected:
	virtual int vertex_setup() = 0;
	virtual int shader_loader(const char *vert_path , const char *frag_path) = 0; 

	const char* API_Name; 
	int ver_major, ver_minor;
	vec2<int> window_size; 
	vec3<int> grid_size; 
	short rendermode; 
	double dt, et; 
};

// OpenGL Render Object 3D -
class renderobject_3D_OGL : public renderobject_3D
{
// FCs
friend class fluidsolver_3; 
public:
	renderobject_3D_OGL(const char *api_name, int v_maj, int v_min, const vec2<int> &w_s, const vec3<int> &g_s, GLFWwindow *winptr);
	~renderobject_3D_OGL();
 
	// Public Overrides of OGL Implementation.
	virtual void shader_pipe(fluidobj_3d *f3obj) override final;
	virtual void render_loop(rend_state dbg) override final;

	// Util 

	enum use_cube
	{
		CUBE_FRONT = 0,
		CUBE_BACK
	};
	
	enum use_fbo
	{
		FBO_CUBE_MS = 0,
		FBO_CUBE_S,
		FBO_DEFAULT
	};

protected:
	// RenderObject Virtual MFunc OVerrides - 
	virtual int vertex_setup() override final;
	virtual int shader_loader(const char *vert_path, const char *frag_path) override final;

	// OGL Specfic MFuncs. 
	void shader_checkCompile(const char *type, int shader);
	void shader_checkLink(int shader_prog); 

	// DVR Setup - 
	void inital_renderstate(); // Called Once Only (in ctor).
	void cube_setup(); 
	void cube_update(); 
	void cube_fbo_setup();
	void cube_fbo_attachAndblit(use_cube tex); // 0 front | 1 back 
	void bindclear_fbo(use_fbo mode); // 0 Cube FBO | 1 Default FBO
	void get_input(const vec2<float> &m);

	// DBG - 
	void print_GL_error(); 
	void calc_FPS(); 
	float get_FPS(short mode);

private:
	// Buffers -
	GLuint Cube_VAO, Cube_Edge_VAO, Quad_VAO;
	GLuint Cube_VBO, Cube_Edge_VBO, Quad_VBO;
	GLuint Quad_EBO, Cube_FBO_s, Cube_FBO_ms, Cube_RBO;

	// Shaders -
	GLuint cube_vert_shader, cube_frag_shader, quad_vert_shader, quad_frag_shader; // 0Cube(vs|fs), 1Quad(vs|fs)
	GLuint cube_shader_prog, quad_shader_prog; // 0Cube, 1Quad
	int cur_shader = 0; // Hacky to avoid changing Base shaderload Parms.  

	// 2D Textures - 
	GLuint tex_CFront_s, tex_CBack_s, tex_ms; 

	// 3D Textures - 
	GLuint tex_dens, tex_vel;

	// Geo Arrays
	GLfloat *cube_vertices, *cube_edge_vertices, *quad_vertices; 
	GLuint *quad_indices; 

	// Cube Transform - 
	matrix_4x4<float> cube_model, cube_view, cube_persp;
	vec3<float> cam_target;
	float zoom;

	// RenderContext (GLFW) Window Pointer
	GLFWwindow *window_ptr = nullptr; 

	// Util Members
	double t0 = 0.0, t1 = 0.0, dt = 0.0, fps = 0.0; 
	uint64_t step = 0;

	// Shader Code Buffers - 
	const char *cube_vert_shader_code, *cube_frag_shader_code, *quad_vert_shader_code, *quad_frag_shader_code;
};



#endif