#ifndef RENDEROBJECT_H
#define RENDEROBJECT_H

#define GLEW_STATIC // Static Linking GLEW and GLFW For now. 
#include <GLEW\glew.h>
#include <GLFW\glfw3.h>

#include "grids2.h"


// Forward Decls - 
class fluidobj_3d; 

// Render State Enum 
enum class rend_state
{
	RENDER_DEBUG,
	RENDER_ACTIVE
};


// Abstract Base Class Interface for API Specfic Render Object. (2D Fluid)
class renderobject_3D
{
// FCs
friend class fluidsolver_3; // Allow Private/Protected Acess of FluidSolver_2. 

public:
	renderobject_2D(const char *api_name, int v_maj, int v_min, int win_x, int win_y, int rmode);
	virtual ~renderobject_2D() = default;

protected:
	virtual int vertex_setup() = 0;
	virtual int shader_loader(const char *vert_path , const char *frag_path) = 0; 
	virtual void shader_pipe(fluidobj_2d *f2obj) = 0;
	virtual void render_loop(rend_state dbg) = 0;

	const char* API_Name; 
	int ver_major, ver_minor;
	int winsize_x, winsize_y; 
	int rendermode; 
	double dt, et; 

};

// OpenGL Render Object 2D -
class renderobject_2D_OGL : public renderobject_2D
{
// FCs
friend class fluidsolver_2; // Allow Private/Protected Acess of FluidSolver_2. 
public:
	renderobject_2D_OGL(const char *api_name, int v_maj, int v_min, int wx, int wy, GLFWwindow *winptr, int rmode);
	~renderobject_2D_OGL();

	// Temp Public Render Getter - 
	void call_ren(rend_state dbg); 
	
protected:
	// RenderObject Virtual MFunc OVerrides - 
	virtual int vertex_setup() override;
	virtual int shader_loader(const char *vert_path, const char *frag_path) override;
	virtual void shader_pipe(fluidobj_2d *f2obj) override;
	virtual void render_loop(rend_state dbg) override;

	// OGL Specfic MFuncs. 
	void shader_checkCompile(const char *type);
	void shader_checkLink(); 

	// DBG - 
	void print_GL_error(); 


private:
	// Buffers -
	GLuint VBO, VAO, EBO;
	// Shaders -
	GLuint vert_shader, frag_shader, shader_prog;
	// Textures - 
	GLuint	tex_dens, tex_vel_u, tex_vel_v, tex_c, tex_vc_u, tex_vc_v, tex_img_rgb;
	// Debug Textures - 
	GLuint tex_preprojvel_u, tex_preprojvel_v, tex_divergence, tex_pressure; 
	// Geo Arrays
	GLfloat *vertices = nullptr; 
	GLuint *indices = nullptr; 

	// RenderContext (GLFW) Window Pointer
	GLFWwindow *window_ptr = nullptr; 

	// Shader Code Buffers - 
	const char *vert_shader_code, *frag_shader_code;
};





#endif