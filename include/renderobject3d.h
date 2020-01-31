#ifndef RENDEROBJECT_H
#define RENDEROBJECT_H

// Prevent Glew and GLFW multiple inclusion. 
#ifndef GLEW_STATIC
#define GLEW_STATIC // Static Linking GLEW and GLFW For now. 
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
	renderobject_3D(const char *api_name, int v_maj, int v_min, const vec2<int> &win_size, const vec3<int> &grd_size, short rmode);
	virtual ~renderobject_3D() = default;

protected:
	virtual int vertex_setup() = 0;
	virtual int shader_loader(const char *vert_path , const char *frag_path) = 0; 
	virtual void shader_pipe(fluidobj_3d *f3obj) = 0;
	virtual void render_loop(rend_state rs) = 0;

	const char* API_Name; 
	int ver_major, ver_minor;
	vec2<int> window_size; 
	vec3<int> grid_size; 
	short rendermode; 
	double dt, et; 

};

// OpenGL Render Object 2D -
class renderobject_3D_OGL : public renderobject_3D
{
// FCs
friend class fluidsolver_3; 
public:
	renderobject_3D_OGL(const char *api_name, int v_maj, int v_min, const vec2<int> &w_s, const vec3<int> &g_s, GLFWwindow *winptr, short rmode);
	~renderobject_3D_OGL();

	// Publicly Callable Start Render. 
	void call_ren(rend_state rs); 
	
protected:
	// RenderObject Virtual MFunc OVerrides - 
	virtual int vertex_setup() override final;
	virtual int shader_loader(const char *vert_path, const char *frag_path) override final;
	virtual void shader_pipe(fluidobj_3d *f3obj) override final;
	virtual void render_loop(rend_state dbg) override final;

	// OGL Specfic MFuncs. 
	void shader_checkCompile(const char *type);
	void shader_checkLink(); 

	// DBG - 
	void print_GL_error(); 


private:
	// Buffers -
	GLuint Quad_VBO, Quad_VAO, Quad_EBO;
	// Shaders -
	GLuint vert_shader, frag_shader, shader_prog;
	// 3D Textures - 
	GLuint tex_dens, tex_vel;
	// Geo Arrays
	GLfloat *vertices = nullptr; 
	GLuint *indices = nullptr; 

	// RenderContext (GLFW) Window Pointer
	GLFWwindow *window_ptr = nullptr; 

	// Shader Code Buffers - 
	const char *vert_shader_code, *frag_shader_code;
};





#endif